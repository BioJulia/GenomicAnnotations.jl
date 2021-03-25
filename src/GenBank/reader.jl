# GenBank Reader

struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
	io::S
end

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
	if input[end-2:end] == ".gz"
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
    feature, posstring = split(strip(line), r" +")
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
Parse footer (sequence) portion of a GenBank file, returning a `String`.
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

    # chromosome = Record{G}()
	chromosome = record

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
		if isheader && occursin(r"FEATURES", line)
            chromosome.header = String(take!(iobuffer))
			isheader = false

        elseif isheader
            linecount == 1 ? print(iobuffer, line) : print(iobuffer, '\n', line)

        # Check if the footer has been reached
		elseif !isheader && !isfooter && (occursin(r"^BASE COUNT", line) || occursin(r"ORIGIN", line))
			# Stop parsing the file when the list of genes is over
            isfooter = true
			iobuffer = IOBuffer()

        ### BODY
		elseif !isheader && !isfooter
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
            if occursin(r"^ {5}\S", line)
                spanning = false
                try
                    feature, locus = parseposition(line)
                catch
                    println(line)
                    @error "parseposition(line) failed at line $linecount"
                end
                addgene!(chromosome, feature, locus)
            elseif !spanning && occursin(r"^ +/", line)
                if occursin(r"=", line)
                    if occursin("=\"", line)
                        (qualifier, content) = match(r"^ +/([^=]+)=\"?([^\"]*)\"?$", line).captures
                        content = String(content)
                    else
                        (qualifier, content) = match(r"^ +/(\S+)=(\S+)$", line).captures
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

					isempty(names(chromosome.genedata)) ?
						(chromosome.genedata[!, Symbol(qualifier)] = Union{Missing, typeof(content)}[content]) :
	                    pushproperty!(chromosome.genes[end], Symbol(qualifier), content)

                else
                    # Qualifiers without a value assigned to them end up here
                    qualifier = split(line, '/')[2]
                    pushproperty!(chromosome.genes[end], Symbol(qualifier), true)
                end
            elseif spanning
                try
                    content = match(r" {21}([^\"]*)\"?$", line)[1]
                catch
                    @warn "Couldn't read content (line $linecount)"
                end
                if line[end] == '"'
                    spanning = false
                end
                if eltype(chromosome.genedata[!, Symbol(qualifier)]).b <: AbstractArray
                    chromosome.genedata[!, Symbol(qualifier)][end][end] = oneline(Base.getproperty(chromosome.genes[end], Symbol(qualifier))[end] * "\n" * content)
                else
                    Base.setproperty!(chromosome.genes[end], Symbol(qualifier), oneline(Base.getproperty(chromosome.genes[end], Symbol(qualifier)) * "\n" * content))
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
    chromosome.name = parseheader(chromosome.header)
    chromosome.sequence = LongDNASeq(filterseq(iobuffer))
    return chromosome
end


"""
    readgbk(input, G::Type = Gene; gunzip = false)

Parse GenBank-formatted file, returning a `Vector{GenBank.Record}`. `input` can be a file path or an `IOStream`. File names ending in ".gz" are assumed to be gzipped and are decompressed. Setting `gunzip` to `true` forces this behaviour.
"""
function readgbk(input)
	collect(open(Reader, input))
end