# GFF3 Reader

Reader = GenomicAnnotations.Reader{GFFFormat, <:BioSequences.Alphabet, <:TranscodingStream}

function Base.iterate(reader::Reader{S}, nextone = Record{Gene,S}()) where S
    if BioGenerics.IO.tryread!(reader, nextone) === nothing
        return nothing
    end
    return nextone, Record{Gene,S}()
end

function BioGenerics.IO.tryread!(reader::Reader, output)
    parsechromosome!(reader.io, output)
end

function parsechromosome!(input, record::Record{<:AbstractGene, S}) where S
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
                TranscodingStreams.unread(input, UInt8.(c for c in "$line\n"))
                return record
            end
            loc = if S === AminoAcidAlphabet
                loc = SpanLocus(parse(Int, sstart):parse(Int, send), ClosedSpan)
            elseif strand == "+"
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
                record.sequence = LongSequence{S}(seq)
            else
                # if line[1] ∉ ['A', 'T', 'G', 'C']
                #     println(line)
                # end
                print(iobuffer, line)
            end
        end
    end
    if !isempty(currentfasta)
        record.sequence = LongSequence{S}(String(take!(iobuffer)))
    end
    record.header = String(take!(header))
    if isempty(record.name)
        return nothing
    else
        return record
    end
end
