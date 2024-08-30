# GenBank Reader

struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    io::S
end

"""
    GenBank.Reader(input::IO)

Create a data reader of the GenBank file format.

```julia
open(GenBank.Reader, "test/example.gbk") do records
    for record in record
        print(record)
    end
end
```
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
Parse footer (sequence) portion of a GenBank file, returning a `String`.
"""
function filterseq(io::IOBuffer)
    line = String(take!(io))
    String(replace(filter(x -> !(isspace(x) || isnumeric(x)), line), r"[^ATGCNatgcn]" => "n"))
end


"""
    parsefeature!(record, geneindex, bytes)

Parse one feature and add it to `record`.
"""
parsefeature!(record, bytes) = parsefeature!(record, lastindex(record.genes), bytes)
function parsefeature!(record, geneindex, bytes)
    newlines = findall(==(UInt8('\n')), bytes)
    fp = 4 + findfirst(==(UInt8(' ')), bytes[6:newlines[1]])
    feature = Symbol(bytes[6:fp])
    loc = Locus(String(bytes[22:(newlines[1]-1)]))
    qualifiers = Dict{Symbol, Any}()
    qualifier = :locus_tag
    addgene!(record, feature, loc, geneindex)
    for p in [(start = x[1], stop = x[2]) for x in zip(newlines[1:end-1], newlines[2:end])]
        line = bytes[p.start+1:p.stop-1]
        if line[22] == UInt8('/')
            qp = findfirst(==(UInt8('=')), line)
            if !isnothing(qp)
                qualifier = Symbol(line[23:qp-1])
                if line[qp+1] == UInt8('"')
                    ### Strings
                    if line[end] != UInt8('"')
                        content = String(line[qp+2:end])
                    else
                        content = String(line[qp+2:end-1])
                    end
                elseif line[qp+1] == UInt8('(')
                    ### Anticodon strings
                    content = String(line[qp+1:end])
                else
                    ### Ints etc.
                    content = String(line[qp+1:end])
                    try
                        tmpcontent = Meta.parse(content)
                        tmpcontent isa Expr && throw(Meta.ParseError)
                        content = tmpcontent
                    catch
                        content = Symbol(content)
                    end
                end
            else
                ### pseudo etc.
                qualifier = Symbol(line[23:end])
                content = true
            end
            ### add feature
            if haskey(qualifiers, qualifier)
                if typeof(qualifiers[qualifier]) isa AbstractVector
                    push!(qualifiers[qualifier], content)
                else
                    qualifiers[qualifier] = [qualifiers[qualifier], content]
                end
            else
                qualifiers[qualifier] = content
            end
        else
            ### Spanning
            if line[end] == UInt8('"')
                content = String(line[22:end-1])
            else
                content = String(line[22:end])
            end
            if haskey(qualifiers, qualifier)
                if qualifiers[qualifier] isa AbstractVector
                    qualifiers[qualifier][end] = oneline(qualifiers[qualifier][end] * "\n" * content)
                else
                    qualifiers[qualifier] = oneline(qualifiers[qualifier] * "\n" * content)
                end
            end
        end
    end
    for (k, v) in qualifiers
        pushproperty!(record, geneindex, k, v)
    end
end


"""
    parsechromosome!(stream::IO, record::Record{G}) where G <: AbstractGene

Parse and return one chromosome entry.
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
    loc = Locus("1..1")


    linecount = 0
    ### HEADER
    while !eof(stream)
        line = readuntil(stream, UInt8('\n'), keep = false)
        linecount += 1
        length(line) <= 1 && continue
        # Catch cases where there's no header
        if linecount == 1 && all(==(UInt8(' ')), line[1:5])
            missingheader = true
            break
        end
        if length(line) >= 8 && line[1:8] == codeunits("FEATURES")
            record.header = String(take!(iobuffer))
            break
        else
            linecount == 1 ? write(iobuffer, line) : write(iobuffer, '\n', line)
        end
    end

    ### BODY
    geneindex = 0
    genebuffer = IOBuffer()
    while !eof(stream)
        line = readuntil(stream, UInt8('\n'), keep = true)
        linecount += 1
        length(line) <= 1 && continue
        if line[6] != UInt8(' ') && all(==(UInt8(' ')), line[1:5])
            if geneindex > 0
                parsefeature!(record, geneindex, take!(genebuffer))
            end
            write(genebuffer, line)
            geneindex += 1
        elseif line[1] != UInt8(' ') && ((length(line) >= 10 && line[1:10] == codeunits("BASE COUNT")) || (length(line) >= 6 && line[1:6] == codeunits("ORIGIN")))
            parsefeature!(record, geneindex, take!(genebuffer))
            break
        elseif line[1] != UInt8(' ')
            continue
        else
            write(genebuffer, line)
        end
    end

    ### FOOTER
    iobuffer = IOBuffer()
    while !eof(stream)
        line = readuntil(stream, UInt8('\n'), keep = true)
        if line == codeunits("//\n")
            break
        end
        write(iobuffer, line)
    end
    record.name = parseheader(record.header)
    record.sequence = LongDNA{4}(filterseq(iobuffer))
    return record
end
