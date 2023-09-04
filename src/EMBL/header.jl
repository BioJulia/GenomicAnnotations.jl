mutable struct HeaderEntry
    # indent::Int
    name
    data
    children::Vector{HeaderEntry}
    HeaderEntry() = new("", "", [])
    HeaderEntry(name, data, children) = new(name, data, children)
end

function Base.empty!(e::HeaderEntry)
    e.name = ""
    e.data = ""
    empty!(e.children)
    nothing
end

function Base.copy(e::HeaderEntry)
    HeaderEntry(e.name, e.data, e.children)
end

mutable struct HeaderLocus
    name
    size
    type
    circular
    division
    date
end
HeaderLocus() = HeaderLocus("", 0, "", false, "", "")
Base.isempty(h::HeaderLocus) = isempty(h.name) && isempty(h.type) && isempty(h.date)

struct Header
    locus::HeaderLocus
    data::Vector{HeaderEntry}
end
Header() = Header(HeaderLocus(), [])
Base.isempty(h::Header) = isempty(h.locus) && isempty(h.data)

struct HeaderEntryView <: AbstractVector{HeaderEntry}
    name
    indices
    parent
end
Base.size(v::HeaderEntryView) = size(v.indices)
Base.getindex(v::HeaderEntryView, I...) = getindex(v.parent[v.indices], I...)
Base.setindex!(v::HeaderEntryView, i, x) = setindex!(view(v.parent, v.indices), i, x)

function parseheaderlocus(s)
    s = s[13:end]
    m = match(r"^(\S+) +(\d+) bp +(...) +(\w+)? +(...) +(\S+)$", s)
    HeaderLocus(m.captures...)
end

indent(line) = findfirst(!=(' '), line) - 1

function parseentry!(current, lines)
    current.name = match(r"\S+", first(lines)[1:12]).match
    current.data = first(lines)[13:end]
    topindent = indent(first(lines))
    length(lines) == 1 && return current
    currentchild = findfirst(line -> line[1:12] != " "^12, lines[2:end])
    if isnothing(currentchild)
        for line in lines[2:end]
            current.data *= "\n" * line[13:end]
        end
        return current
    else
        currentchild += 1
    end
    for line in lines[2:currentchild-1]
        current.data *= "\n" * line[13:end]
    end
    # current.data = GenomicAnnotations.oneline(current.data)
    currentchildindent = indent(lines[currentchild])
    nextchild = findfirst(line -> indent(line) == currentchildindent, view(lines, currentchild+1:lastindex(lines)))
    if isnothing(nextchild)
        nextchild = lastindex(lines)+1
    else
        nextchild += currentchild + 1
    end
    while currentchild <= lastindex(lines)
        currentchildentry = HeaderEntry()
        if length(currentchild:nextchild-1) < 1
            return current
        end
        parseentry!(currentchildentry, view(lines, currentchild:nextchild-1))
        push!(current.children, copy(currentchildentry))
        currentchild = nextchild
        nextchild = findfirst(line -> indent(line) == currentchildindent, view(lines, currentchild+1:lastindex(lines)))
        if isnothing(nextchild)
            break
        else
            nextchild = nextchild + 1
        end
    end
    return current
end

function parseheader(::Type{Header}, s)
    lines = split(s, '\n')
    entries = HeaderEntry[]
    ### LOCUS
    headerlocus = parseheaderlocus(lines[1])
    ### Other entries
    i = 2
    while i <= lastindex(lines)
        current = HeaderEntry()
        nextmain = findfirst(s -> s[1] != ' ', lines[i+1:end])
        if isnothing(nextmain)
            parseentry!(current, view(lines, i:lastindex(lines)))
            push!(entries, copy(current))
            i = lastindex(lines) + 1
            break
        else
            parseentry!(current, view(lines, i:i+nextmain-1))
            push!(entries, copy(current))
            # empty!(current)
            i = i + nextmain
        end
    end
    Header(headerlocus, entries)
end

parseheader(::Type{String}, s) = s

function Base.show(io::IO, h::HeaderLocus)
    println(io, "LOCUS       ",
        rpad(h.name, 16, " "),
        lpad(h.size, 12, " "),
        " bp",
        lpad(h.type, 6, " "),
        lpad(h.circular, 12, " "),
        lpad(h.division, 5, " "),
        " ",
        h.date)
end

Base.show(io::IO, mime::MIME"text/plain", A::AbstractArray{HeaderEntry}) = show(io, A)
function Base.show(io::IO, A::AbstractArray{HeaderEntry})
    for he in A
        showindent(io, he, 1)
    end
end
Base.show(io::IO, he::HeaderEntry) = showindent(io, he, 1)
function showindent(io::IO, he::HeaderEntry, level::Int)
    indent = level == 1 ? 0 : level
    print(io, rpad(" "^indent * he.name, 12, " "))
    data = split(he.data, "\n")
    println(io, data[1])
    for d in data[2:end]
        println(io, " "^12, d)
    end
    for child in he.children
        showindent(io, child, level+1)
    end
end

Base.show(io::IO, mime::MIME"text/plain", header::Header) = show(io, header)
function Base.show(io::IO, header::Header)
    ### LOCUS
    show(io, header.locus)
    ### Other entries
    for entry in header.data
        showindent(io, entry, 1)
    end
end

Base.getindex(h::Header, s::Symbol) = getindex(h, string(s))
function Base.getindex(h::Header, s::String)
    if s == "LOCUS"
        return h.locus
    end
    indices = findall(e -> e.name == s, h.data)
    HeaderEntryView(s, indices, h.data)
end

Base.getindex(h::HeaderEntry, s::Symbol) = getindex(h, string(s))
function Base.getindex(h::HeaderEntry, s::String)
    indices = findall(e -> e.name == s, h.children)
    HeaderEntryView(s, indices, h.children)
end

function Base.setindex!(A::HeaderEntryView, v::AbstractString, i::Int)
    if i in eachindex(A)
        A[i].data = v
    elseif i == lastindex(A) + 1
        push!(A.parent, HeaderEntry(A.name, v, []))
    end
end

iscircular!(h::Header, b::Bool = true) = h.locus.circular = b ? "circular" : "linear"

parsename(header::GenBank.Header) = header.locus.name