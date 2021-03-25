# GFF3 Reader

struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    io::S
end

"""
    GFF.Reader(input::IO)

Create a data reader of the GFF3 file format.
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
    parsechromosome!(reader.io, output)
end

function parsechromosome!(input, record::Record{G}) where G <: AbstractGene
    record.genedata[!, :source] = Union{Missing, String}[]
    record.genedata[!, :score] = Union{Missing, Float64}[]
    record.genedata[!, :phase] = Union{Missing, Int}[]
    iobuffer = IOBuffer()
    isheader = true
    isfooter = false
    qualifier = String("")
    content = String("")
    header = IOBuffer()
    currentfasta = ""
    linecount = 0
    while !eof(input)
        linecount += 1
        line = readline(input)
        if length(line) == 0
            continue
        end
        # Store the header
        if line[1] == '#' && line != "##FASTA"
            isheader = true
            println(header, line)
        else
            isheader = false
        end
        if line == "##FASTA"
            isfooter = true
        elseif !isfooter && !isheader
            (seqid, source, feature, sstart, send, score, strand, phase, attributes) = split(line, '\t')
            if isempty(record.name)
                record.header = String(take!(header))
                record.name = seqid
            elseif seqid != record.name
                ### The following contig lacks a header
                TranscodingStreams.unread(input, UInt8.(c for c in line))
                return record
            end
            locus = Locus(parse(Int, sstart):parse(Int, send), strand[1])
            addgene!(record, Symbol(feature), locus)
            if source != "." && ismissing(record.genes[end].source)
                record.genes[end].source = source
            end
            if score != "." && ismissing(record.genes[end].score)
                record.genes[end].score = parse(Float64, score)
            end
            if phase != "."
                record.genes[end].phase = parse(Int, phase)
            end
            if attributes != "."
                for attribute in split(attributes, ';')
                    isempty(attribute) && continue
                    qualifier, values = split(attribute, '=')
                    values = split(values, ',')
                    for value in values
                        if occursin(r"^\d+$", value)
                            value = parse(Int, value)
                        else
                            value = String(value)
                        end
                        pushproperty!(record.genes[end], Symbol(qualifier), value)
                    end
                end
            end
        elseif isfooter
            if line[1] == '>'
                currentfasta = line[2:end]
                seq = String(take!(iobuffer))
                record.sequence = LongDNASeq(seq)
            else
                if line[1] ∉ ['A', 'T', 'G', 'C']
                    println(line)
                end
                print(iobuffer, line)
            end
        end
    end
    if !isempty(currentfasta)
        record.sequence = LongDNASeq(String(take!(iobuffer)))
    end
    record.header = headerstring = String(take!(header))
    if isempty(record.name)
        return nothing
    else
        return record
    end
end
