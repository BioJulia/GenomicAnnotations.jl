## GTF/GFF2 Writer

struct Writer{S <: TranscodingStream} <: BioGenerics.IO.AbstractWriter
    output::S
end

function BioGenerics.IO.stream(writer::Writer)
    return writer.output
end

"""
    GTF.Writer(output::IO; width=70)

Create a data writer of the GFF file format.
```julia
open(GTF.Writer, outfile) do writer
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
    printgtf(writer.output, record)
end

function Base.write(writer::Writer, records::AbstractVector{Record})
    for record in records
        printgtf(writer.output, record)
    end
end


function gtfstring(gene::Gene)
    buf = IOBuffer()
    firstattribute = true
    for field in keys(parent(gene).genedata[index(gene), :])
        field in [:source, :score, :phase] && continue
        v = parent(gene).genedata[index(gene), field]
        if !ismissing(v)
            if firstattribute
                firstattribute = false
            else
                print(buf, "; ")
            end
            if v isa AbstractVector
                for i in eachindex(v)
                    i > 1 && print(buf, "; ")
                    print(buf, field, " \"", oneline(v[i]), "\"")
                end
            else
                print(buf, field, " \"", oneline(v), "\"")
            end
        end
    end
    if iscompound(gene)
        s = String(take!(buf))
        res = IOBuffer()
        for loc in locus(gene)
            println(res, join([parent(gene).name,
                get(gene, :source, "."),
                feature(gene),
                loc.start,
                loc.stop,
                get(gene, :score, "."),
                loc.strand,
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
    printgtf(io::IO, chr)
    printgtf(path::AbstractString, chr)

Print `chr` in GTF/GFF2 format.
"""
printgtf(filepath::AbstractString, chrs) = printgtf(open(filepath, "w"), chrs)
printgtf(io::IO, chr::Record) = printgtf(io, [chr])
function printgtf(io::IO, chrs::AbstractVector{Record{G}}) where G <: AbstractGene
    iobuffer = IOBuffer()
    ### GTF files do not have a header
    ### Body
    for chr in chrs
        for gene in chr.genes
            print(iobuffer, gtfstring(gene))
        end
    end
    ### GTF files do not have a footer
    print(io, String(take!(iobuffer)))
end
