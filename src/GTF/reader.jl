# GTF/GFF2 Reader

Reader = GenomicAnnotations.Reader{GTFFormat, <:BioSequences.Alphabet, <:TranscodingStream}

function Base.iterate(reader::Reader{S}, nextone = Record{Gene,S}()) where S
    if BioGenerics.IO.tryread!(reader, nextone) === nothing
        return nothing
    end
    return nextone, Record{Gene,S}()
end

function BioGenerics.IO.tryread!(reader::Reader{S}, output) where S
    parsechromosome!(reader.io, output)
end

function parsechromosome!(input, record::Record)
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
        if line[1] == '#' && line != "##FASTA"
            isheader = true
            println(header, line)
            continue
        else
            isheader = false
        end
        (seqid, source, feature, sstart, send, score, strand, phase, attributes) = split(line, '\t')
        if isempty(record.name)
            record.header = String(take!(header))
            record.name = seqid
        elseif seqid != record.name
            ### The following contig lacks a header
            TranscodingStreams.unread(input, UInt8.(c for c in "$line\n"))
            return record
        end
        loc = if strand == "+"
            SpanLocus(parse(Int, sstart):parse(Int, send), ClosedSpan)
        else
            Complement(SpanLocus(parse(Int, sstart):parse(Int, send), ClosedSpan))
        end
        addgene!(record, Symbol(feature), loc)
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
            for attribute in split(attributes, r" ?; ?")
                isempty(attribute) && continue
                space_pos = findfirst(==(' '), attribute)
                qualifier = attribute[1:(space_pos-1)]
                value = attribute[(space_pos+1):end]
                if occursin(r"^\d+$", value)
                    value = parse(Int, value)
                elseif occursin(r"^\".+\"$", value)
                    value = String(value[2:end-1])
                else
                    value = String(value)
                end
                pushproperty!(record.genes[end], Symbol(qualifier), value)
            end
        end
    end
    record.header = String(take!(header))
    if isempty(record.name)
        return nothing
    else
        return record
    end
end
