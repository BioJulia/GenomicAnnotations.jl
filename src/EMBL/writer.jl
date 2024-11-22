# GenBank Writer

struct Writer{S <: TranscodingStream} <: BioGenerics.IO.AbstractWriter
    output::S
end

function BioGenerics.IO.stream(writer::Writer)
    return writer.output
end

"""
    GenBank.Writer(output::IO; width=70)

Create a data writer of the GenBank file format.
```julia
open(GenBank.Writer, outfile) do writer
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

function Base.write(writer::Writer, record::Record, header = Dict())
    printembl(writer.output, record, header)
end

function Base.write(writer::Writer, records::AbstractVector{Record}, header = Dict())
    for record in records
        printembl(writer.output, record, header)
    end
end

"""
    printgbk(io::IO, chr)
    printgbk(path::AbstractString, chr)

Print `chr` in GenBank format.
"""
function printembl(path::AbstractString, chrs, header = Dict())
    open(path, "w") do f
        printembl(f, chrs, header)
    end
end
function printembl(chrs::AbstractVector{C}, header = Dict()) where {C <: Record}
    io = IOBuffer()
    printembl(io, chrs, header)
end
function printembl(io::IO, chrs::AbstractVector{C}, header = Dict()) where {C <: Record}
    for chr in chrs
        printembl(io, chr, header)
    end
    return io
end
function printembl(chr::C, header = Dict()) where {C <: Record}
    io = IOBuffer()
    printembl(io, chr, header)
end
function printembl(io::IO, chr::C, header = Dict()) where {C <: Record}
    printheader(io, chr, header)
    ### print features
    for gene in chr.genes
        show_embl_gene(io, gene)
    end
    # println(io, chr.genes)
    println(io, "XX")
    println(io, "SQ   Sequence $(length(chr.sequence)) BP; $(count(==(DNA_A), chr.sequence)) A; $(count(==(DNA_C), chr.sequence)) C; $(count(==(DNA_G), chr.sequence)) G; $(count(==(DNA_T), chr.sequence)) T; 0 other;")
    formatsequence(chr.sequence, io)
    println(io, "//")
    return io
end


function appendstring_embl(field, v::Bool)
    return "\n" * rpad("FT", 21, ' ') * "/$field"
end
function appendstring_embl(field, v::Union{Number, Symbol})
    return "\n" * rpad("FT", 21, ' ') * "/$field=" * string(v)
end
function appendstring_embl(field, v)
    v = string("\"", v, "\"")
    v = multiline(v, field; margin = 21)
    v = replace(v, "\n" => "\n" * rpad("FT", 21, ' '))
    return "\n" * rpad("FT", 21, ' ') * "/$field=" * v
end

function show_embl_gene(io::IO, gene::AbstractGene)
    buf = IOBuffer()
    print(buf, "FT   " * rpad(string(feature(gene)), 16, ' '))
    print(buf, locus(gene))
    for field in names(parent(gene).genedata)
        field in [:feature, :locus] && continue
        v = parent(gene).genedata[index(gene), field]
        if !ismissing(v)
            if v isa AbstractVector
                for i in eachindex(v)
                    print(buf, appendstring_embl(field, v[i]))
                end
            else
                print(buf, appendstring_embl(field, v))
            end
        end
    end
    println(io, String(take!(buf)))
end


function printheader(io, chr, header = Dict())
    println(io, "ID   XXX; XXX; circular; XXX; XXX; XXX; XXX.")
    println(io, "XX")
    for k in [:project, :organism_species, :organism_classification, :description, :accession, :name]
        !haskey(header, k) && continue
        v = get(header, k, "")
        if k == :project
            println(io, "PR   Project:", v)
        elseif k == :date
            println(io, "DT   ", v[1].date, " (Rel. ", v[1].release, ", Created")
            println(io, "DT   ", v[2].date, " (Rel. ", v[2].release, ", Last updated, Version ", v[2].version, ")")
        elseif k == :organism_species
            println(io, "OS   ", v)
        elseif k == :organism_classification
            if length(v) > 80
                w = String[]
                I = findall(==(' '), v)
                i = I[findlast(<=(80), I)]
                push!(w, v[1:i-1])
                push!(w, v[i+1:end])
                foreach(v -> println(io, "OC   ", v), w)
            else
                println(io, "OC   ", v)
            end
        elseif k == :description
            println(io, "DE   ", v)
        elseif k == :accession
            println(io, "AC   ", v)
        elseif k == :name
            println(io, "AC * _", v)
        end
        println(io, "XX")
    end
    println(io, "FH   Key           Location/Qualifiers")
    println(io, "FH")
end


function formatsequence(sequence, io = IOBuffer)
    p = length(string(length(sequence))) + 2
    if length(sequence) > 60
        intervals = [i:i+60 for i in range(1; step = 60, stop = length(sequence)-60)]
        for interval in intervals
            println(io, "     ", sequence[interval[1:10]],
                " ", sequence[interval[11:20]], " ", sequence[interval[21:30]],
                " ", sequence[interval[31:40]], " ", sequence[interval[41:50]],
                " ", sequence[interval[51:60]], lpad(string(last(interval)), 10, ' '))
        end
    else
        intervals = [1:1]
    end
    i = intervals[end].stop
    if i <= length(sequence)
        print(io, "     ")
        j = 0
        while i+j <= length(sequence)
            print(io, sequence[i+j])
            (j+1) % 10 == 0 && print(io, " ")
            j += 1
        end
    end
    println(io)
    return io
end
