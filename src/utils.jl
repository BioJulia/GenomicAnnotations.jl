"""
    oneline(s)

Remove newlines, replacing them with spaces when it seems appropriate.
"""
function oneline(v::AbstractString)
    R = split(v, '\n')
    buf = IOBuffer()
    for (i, r) in enumerate(R)
        if i != lastindex(R) && (r[end] == '-' || occursin(' ', r))
            print(buf, strip(r), " ")
        else
            print(buf, r)
        end
    end
    String(take!(buf))
end
oneline(v::Any) = v


function findbreak(v, t; width = 80, margin = 21)
    i = t == width - margin + 1 ?
        findlast('\n', v) :
        1
    j = i+t-1
    j >= lastindex(v) && return nothing
    x = findlast(c -> c == ' ' || c == '-', v[j-t+1:j])
    if isnothing(x)
        return j:j+1
    elseif v[x] == ' '
        x += j - t
        return x-1:x+1
    else
        x += j - t
        return x:x+1
    end
end


"""
    multiline(v, s::Symbol; width = 80, margin = 21, first_line_margin = 4)

Return a `String` with newlines and spaces added so that it conforms to the specified line width, with the margin subtracted. For example, the line width for the GenBank format is 80, with a margin of 21, whereas for the EMBL format the margin is 19.
"""
function multiline(v, s; width = 80, margin = 21, first_line_margin = 3)
    v = string(v)
    s = string(s)
    if length(v) + length(s) > width - margin - first_line_margin + 1 && !occursin('\n', v)
        b = findbreak(v, width - margin - first_line_margin + 1 - length(s))
        v = v[1:b.start] * "\n" * v[b.stop:end]
        b = findbreak(v, width - margin + 1; width, margin)
        while !isnothing(b) && b.start < length(v)
            v = v[1:b.start] * "\n" * v[b.stop:end]
            b = findbreak(v, width - margin + 1; width, margin)
        end
    end
    return v
end


"""
    readgbk(input)

Parse GenBank-formatted file, returning a `Vector{Record}`. File names ending in ".gz" are assumed to be gzipped and are decompressed.
"""
function readgbk(input)
    collect(open(GenBank.Reader, input))
end


"""
    readembl(input)

Parse EMBL-formatted file, returning a `Vector{Record}`. File names ending in ".gz" are assumed to be gzipped and are decompressed.
"""
function readembl(input)
    collect(open(EMBL.Reader, input))
end


"""
    readgff(input)

Parse GFF3-formatted file, returning a `Vector{Record}`. File names ending in ".gz" are assumed to be gzipped and are decompressed.
"""
function readgff(input)
    collect(open(GFF.Reader, input))
end


"""
    readgff(input)

Parse GTF/GFF2-formatted file, returning a `Vector{Record}`. File names ending in ".gz" are assumed to be gzipped and are decompressed.
"""
function readgtf(input)
    collect(open(GTF.Reader, input))
end


"""
    relative_position(gene::AbstractGene, point::Symbol = :start)

Return the relative position of `gene` on `parent(gene)`, on the interval (0,1]. The keyword `point` determines which of the `:start`, `:middle`, or `:stop` should be counted.
"""
relative_position(gene::AbstractGene, point = :start) = relative_position(parent(gene).sequence, locus(gene), point; iscomplement = iscomplement(gene))
function relative_position(chrseq, loc::AbstractLocus, point = :start; iscomplement = false)
    if (point == :start && !iscomplement) || (point == :stop && iscomplement)
        return loc.start / length(chrseq)
    elseif (point == :stop && !iscomplement) || (point == :start && iscomplement)
        return loc.stop / length(chrseq)
    elseif point == :middle
        p1 = loc.start / length(chrseq)
        p2 = loc.stop / length(chrseq)
        return mod1(atan((sinpi(2p1) + sinpi(2p2)) / 2, (cospi(2p1) + cospi(2p2)) / 2) / 2π, 1)
    else
        error("`point` must be one of `:start`, `:middle`, or `:stop`")
    end
end