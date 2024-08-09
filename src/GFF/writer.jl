## GFF3 Writer

struct Writer{S <: TranscodingStream} <: BioGenerics.IO.AbstractWriter
    output::S
end

function BioGenerics.IO.stream(writer::Writer)
    return writer.output
end

"""
    GFF.Writer(output::IO; width=70)

Create a data writer of the GFF file format.
```julia
open(GFF.Writer, outfile) do writer
    write(writer, genome)
end
```
"""
function Writer(output::IO)
    if output isa TranscodingStream
        return Writer{typeof(output)}(stream)
    else
        stream = TranscodingStreams.NoopStream(output)
        return Writer{typeof(stream)}(stream)
    end
end

function Base.write(writer::Writer, record::Record)
    printgff(writer.output, record)
end


function gffstring(gene::Gene)
    buf = IOBuffer()
    firstattribute = true
    for field in names(parent(gene).genedata)
        field in [:source, :score, :phase] && continue
        v = parent(gene).genedata[index(gene), field]
        if !ismissing(v)
            if firstattribute
                firstattribute = false
            else
                print(buf, ";")
            end
            if v isa AbstractVector
                print(buf, field, "=")
                for i in eachindex(v)
                    print(buf, oneline(v[i]))
                    i == lastindex(v) ? print(buf, ";") : print(buf, ",")
                end
            else
                print(buf, field, "=", oneline(v))
            end
        end
    end
    if locus(gene).position isa Vector
        s = String(take!(buf))
        res = IOBuffer()
        for pos in locus(gene).position
            println(res, join([parent(gene).name,
                get(gene, :source, "."),
                feature(gene),
                pos.start,
                pos.stop,
                get(gene, :score, "."),
                locus(gene).strand,
                get(gene, :phase, "."),
                s], '\t'))
        end
        String(take!(res))
    else
        join([parent(gene).name,
            get(gene, :source, "."),
            feature(gene),
            locus(gene).start,
            locus(gene).stop,
            get(gene, :score, "."),
            locus(gene).strand,
            get(gene, :phase, "."),
            String(take!(buf))], '\t') * "\n"
    end
end

"""
    printgff(io::IO, chr)
    printgff(path::AbstractString, chr)

Print `chr` in GFF3 format.
"""
printgff(filepath::AbstractString, chrs) = printgff(open(filepath, "w"), chrs)
printgff(io::IO, chr::Record) = printgff(io, [chr])
function printgff(io::IO, chrs::AbstractVector{Record{G}}) where G <: AbstractGene
    iobuffer = IOBuffer()
    ### Header
    if chrs[1].header[1:15] == "#gff-version 3\n"
        print(iobuffer, chrs[1].header)
    else
        println(iobuffer, "##gff-version 3")
    end
    ### Body
    for chr in chrs
        for gene in chr.genes
            print(iobuffer, gffstring(gene))
        end
    end
    ### Footer
    if !all(isempty(chr.sequence) for chr in chrs)
        println(iobuffer, "##FASTA")
        for chr in chrs
            println(iobuffer, ">", chr.name)
            for s in Iterators.partition(chr.sequence, 80)
                println(iobuffer, join(s))
            end
        end
    end
    print(io, String(take!(iobuffer)))
end
