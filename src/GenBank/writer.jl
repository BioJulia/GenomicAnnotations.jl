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

function Base.write(writer::Writer, record::Record)
    printgbk(writer.output, record)
end

"""
    printgbk(io::IO, chr)
    printgbk(path::AbstractString, chr)

Print `chr` in GenBank format.
"""
function printgbk(path::AbstractString, chrs)
    open(path, "w") do f
        printgbk(f, chrs)
    end
end
function printgbk(chrs::AbstractVector{C}) where {C <: Record}
    io = IOBuffer()
    printgbk(io, chrs)
end
function printgbk(io::IO, chrs::AbstractVector{C}) where {C <: Record}
    for chr in chrs
        printgbk(io, chr)
    end
    return io
end
function printgbk(chr::C) where {C <: Record}
    io = IOBuffer()
    printgbk(io, chr)
end
function printgbk(io::IO, chr::C) where {C <: Record}
    if !isempty(chr.header) && !occursin(r"^#", chr.header)
        println(io, chr.header)
    else
        ### GFF3 header
        println(io, "LOCUS       unknown", lpad(length(chr.sequence), 10, ' '), " bp     DNA")
    end
    println(io, rpad("FEATURES", 21, ' '), "Location/Qualifiers")
    println(io, chr.genes)
    println(io, "ORIGIN")
    formatsequence(chr.sequence, io)
    println(io)
    println(io, "//")
    return io
end


function formatsequence(sequence, io = IOBuffer)
    p = length(string(length(sequence))) + 2
    if length(sequence) > 60
        intervals = [i:i+60 for i in range(1; step = 60, stop = length(sequence)-60)]
        for interval in intervals
            println(io, lpad(string(first(interval)), p, ' '), " ", sequence[interval[1:10]],
                " ", sequence[interval[11:20]], " ", sequence[interval[21:30]],
                " ", sequence[interval[31:40]], " ", sequence[interval[41:50]],
                " ", sequence[interval[51:60]])
        end
    else
        intervals = [1:1]
    end
    i = intervals[end].stop
    if i <= length(sequence)
        print(io, lpad(i, p, ' '), " ")
        j = 0
        while i+j <= length(sequence)
            print(io, sequence[i+j])
            (j+1) % 10 == 0 && print(io, " ")
            j += 1
        end
    end
    return io
end
