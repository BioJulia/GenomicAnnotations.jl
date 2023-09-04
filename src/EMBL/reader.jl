# EMBL Reader

struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    io::S
end

"""
    EMBL.Reader(input::IO)

Create a data reader of the EMBL file format.
"""
function Reader(input::IO)
    if input isa TranscodingStream
        Reader(input)
    else
        stream = TranscodingStreams.NoopStream(input)
        Reader(stream)
    end
end

function Base.eltype(::Type{<:Reader})
    return Record
end

function Base.close(reader::Reader)
    close(reader.io)
end

function Base.open(::Type{<:Reader}, input::AbstractString)
    if input[end-2:end] == ".gz"
        return Reader(GzipDecompressorStream(open(input)))
    else
        return Reader(TranscodingStreams.NoopStream(open(input)))
    end
end

function BioGenerics.IO.stream(reader::Reader)
    reader.io
end

function Base.iterate(reader::Reader, nextone::Record = Record())
    if BioGenerics.IO.tryread!(reader, nextone) === nothing
        return nothing
    end
    return nextone, Record()
end

function BioGenerics.IO.tryread!(reader::Reader, output)
    record = parsechromosome!(reader.io, output)
    return record
end

"""
Return the LOCUS entry of the header.
"""
function parseheader(header::String)
    lines = split(header, "\n")
    firstline = lines[findfirst(line -> occursin("LOCUS", line), lines)]
    locus = split(firstline, r" +")[2]
    return locus
end


"""
Parse lines encoding genomic position, returning the feature as a `Symbol`, and an instance of `Locus`.
"""
function parseposition(line::String)
    FT, feature, posstring = split(strip(line), r" +")
    if occursin(r"(\.\.|\^)", posstring)
        position = UnitRange(parse.(Int, filter.(c->isnumeric(c), split(posstring, r"(\.\.|\^)(.*(\.\.|\^))?")))...)
        strand = occursin("complement", posstring) ? '-' : '+'
    else
        position = parse(Int, posstring)
        position = position:position
        strand = '.'
    end
    complete_left = !occursin('<', posstring)
    complete_right = !occursin('>', posstring)
    order = Vector{UnitRange{Int}}()
    join = occursin("join", posstring)
    if join || occursin("order", posstring)
        for m in eachmatch(r"\d+(\.\.|\^)\d+", posstring)
            r = Meta.parse.(split(m.match, r"(\.\.|\^)"))
            push!(order, r[1]:r[2])
        end
    end
    return Symbol(feature), Locus(position, strand, complete_left, complete_right, order, join)
end


"""
Parse footer (sequence) portion of a EMBL file, returning a `String`.
"""
function filterseq(io::IOBuffer)
    line = String(take!(io))
    String(replace(filter(x -> !(isspace(x) || isnumeric(x)), line), r"[^ATGCNatgcn]" => "n"))
end


"""
    parsechromosome(stream::IO)

Parse and return one chromosome entry, and the line number that it ends at.
"""
function parsechromosome!(stream::IO, record::Record{G}) where G <: AbstractGene
    eof(stream) && return nothing
    iobuffer = IOBuffer()
    isheader = true
    isfooter = false
    spanning = false
    position_spanning = false
    qualifier = String("")
    content = String("")
    header = ""

    feature = :source
    locus = Locus()


    linecount = 0
    while !eof(stream)
        line = readline(stream)
        linecount += 1

        if length(line) == 0
            continue
        end

        # Catch cases where there's no header
        if linecount == 1 && occursin(r"gene", line)
            isheader = false
        end

        ### HEADER
        if isheader && occursin(r"^FT", line)
            record.header = String(take!(iobuffer))
            isheader = false

        elseif isheader
            linecount == 1 ? print(iobuffer, line) : print(iobuffer, '\n', line)

        # Check if the footer has been reached
        elseif !isheader && !isfooter && (occursin(r"^SQ", line))
            # Stop parsing the file when the list of genes is over
            isfooter = true
            iobuffer = IOBuffer()
        end

        ### BODY
        if !isheader && !isfooter
            if position_spanning && occursin(r"  /", line)
                position_spanning = false
                spanningline = filter(x -> x != ' ', String(take!(iobuffer)))
                try
                    feature, locus = parseposition(spanningline)
                catch
                    println(spanningline)
                    println(line)
                    @error "parseposition(spanningline) failed at line $linecount"
                end
            elseif position_spanning
                print(iobuffer, line)
            end
            if occursin(r"^FT   \S", line)
                spanning = false
                try
                    feature, locus = parseposition(line)
                catch
                    println(line)
                    @error "parseposition(line) failed at line $linecount"
                end
                addgene!(record, feature, locus)
            elseif !spanning && occursin(r"^FT +/", line)
                if occursin(r"=", line)
                    if occursin("=\"", line)
                        (qualifier, content) = match(r"^FT +/([^=]+)=\"?([^\"]*)\"?$", line).captures
                        content = String(content)
                    else
                        (qualifier, content) = match(r"^FT +/(\S+)=(\S+)$", line).captures
                        try
                            tmpcontent = Meta.parse(content)
                            tmpcontent isa Expr && throw(Meta.ParseError)
                            content = tmpcontent
                        catch
                            content = Symbol(content)
                        end
                    end

                    if occursin(r"\".*[^\"]$", line)
                        spanning = true
                    end

                    isempty(names(record.genedata)) ?
                        (record.genedata[!, Symbol(qualifier)] = Union{Missing, typeof(content)}[content]) :
                        pushproperty!(record.genes[end], Symbol(qualifier), content)

                else
                    # Qualifiers without a value assigned to them end up here
                    qualifier = split(line, '/')[2]
                    pushproperty!(record.genes[end], Symbol(qualifier), true)
                end
            elseif spanning
                try
                    content = match(r"FT {19}([^\"]*)\"?$", line)[1]
                catch
                    @warn "Couldn't read content (line $linecount)"
                end
                if line[end] == '"'
                    spanning = false
                end
                if eltype(record.genedata[!, Symbol(qualifier)]).b <: AbstractArray
                    record.genedata[!, Symbol(qualifier)][end][end] = oneline(Base.getproperty(record.genes[end], Symbol(qualifier))[end] * "\n" * content)
                else
                    Base.setproperty!(record.genes[end], Symbol(qualifier), oneline(Base.getproperty(record.genes[end], Symbol(qualifier)) * "\n" * content))
                end
            end

        ### FOOTER
        elseif isfooter
            if line == "//"
                break
            end
            if occursin(r"^ ", line)
                print(iobuffer, line)
            end
        end
    end
    if isempty(record.header) && isempty(record.genes) && isempty(record.sequence)
        return nothing
    end
    record.name = "" #parseheader(record.header)
    record.sequence = LongDNA{4}(filterseq(iobuffer))
    return record
end
