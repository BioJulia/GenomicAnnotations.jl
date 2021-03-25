# GFF3 Reader

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
    parsechromosome!(reader.io, output)
end

function parsechromosome!(input, output::Record{G}) where G <: AbstractGene
    chr = output
    chr.genedata[!, :source] = Union{Missing, String}[]
    chr.genedata[!, :score] = Union{Missing, Float64}[]
    chr.genedata[!, :phase] = Union{Missing, Int}[]
    iobuffer = IOBuffer()
    isheader = true
    isfooter = false
    qualifier = String("")
    content = String("")
    header = IOBuffer()
    currentchr = 0
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
            # println("", line)
            (seqid, source, feature, sstart, send, score, strand, phase, attributes) = split(line, '    ')
            if isempty(chr.name)
                chr.header = String(take!(header))
                chr.name = seqid
            elseif seqid != chr.name
                ### The following contig lacks a header
                TranscodingStreams.unread(input, UInt8.(c for c in line))
                return chr
            end
            locus = Locus(parse(Int, sstart):parse(Int, send), strand[1])
            addgene!(chr, Symbol(feature), locus)
            if source != "." && ismissing(chr.genes[end].source)
                chr.genes[end].source = source
            end
            if score != "." && ismissing(chr.genes[end].score)
                chr.genes[end].score = parse(Float64, score)
            end
            if phase != "."
                chr.genes[end].phase = parse(Int, phase)
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
                        pushproperty!(chr.genes[end], Symbol(qualifier), value)
                    end
                end
            end
        elseif isfooter
            if line[1] == '>'
                currentfasta = line[2:end]
                seq = String(take!(iobuffer))
                chr.sequence = LongDNASeq(seq)
            else
                if line[1] ∉ ['A', 'T', 'G', 'C']
                    println(line)
                end
                print(iobuffer, line)
            end
        end
    end
    if !isempty(currentfasta)
        chr.sequence = LongDNASeq(String(take!(iobuffer)))
    end
    chr.header = headerstring = String(take!(header))
    if isempty(chr.name)
        return nothing
    else
        return chr
    end
end

"""
    readgff(input, G::Type = Gene; gunzip = false)

Parse GFF3-formatted file, returning a `Vector{Record}`. `input` can be a file path or an `IOStream`. File names ending in ".gz" are assumed to be gzipped and are decompressed. Setting `gunzip` to `true` forces this behaviour.
The type of `AbstractGene` to be used can be specified with `G`, though currently the only option is `Gene`.
"""
function readgff(input)
    collect(open(input))
end
