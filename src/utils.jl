"""
    oneline(s)

Remove newlines, replacing them with spaces when it seems appropriate.
"""
function oneline(v::AbstractString)
    R = split(v, '\n')
    buf = IOBuffer()
    for (i, r) in enumerate(R)
        if i != lastindex(R) && (r[end] == '-' || occursin(' ', r))
            print(buf, r, " ")
        else
            print(buf, r)
        end
    end
    String(take!(buf))
end
oneline(v::Any) = v


function findbreak(v, t)
    i = t == 59 ?
        findlast('\n', v) :
        1
    j = i+t-1
    j >= lastindex(v) && return nothing
    x = findlast(c -> c == ' ' || c == '-', v[1:j])
    if isnothing(x)
        return j:j+1
    elseif x == ' '
        return x-1:x+1
    else
        return x:x+1
    end
end


"""
    multiline(v, s::Symbol)

Return a `String` with newlines and spaces added so that it conforms to the GenBank line width.
"""
function multiline(v, s)
    v = string(v)
    s = string(s)
    if length(v) + length(s) > 55 && !occursin('\n', v)
        b = findbreak(v, 55 - length(s))
        v = v[1:b.start] * "\n" * v[b.stop:end]
        b = findbreak(v, 59)
        while !isnothing(b) && b.start < length(v)
            v = v[1:b.start] * "\n" * v[b.stop:end]
            b = findbreak(v, 59)
        end
    end
    return v
end


"""
    readgbk(input)

Parse GenBank-formatted file, returning a `Vector{GenBank.Record}`. File names ending in ".gz" are assumed to be gzipped and are decompressed.
"""
function readgbk(input)
    collect(open(GenBank.Reader, input))
end


"""
    readgff(input)

Parse GFF3-formatted file, returning a `Vector{Record}`. File names ending in ".gz" are assumed to be gzipped and are decompressed.
"""
function readgff(input)
    collect(open(GFF.Reader, input))
end
