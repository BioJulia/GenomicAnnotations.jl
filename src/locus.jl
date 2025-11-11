abstract type AbstractLocus end

abstract type AbstractDescriptor end
"""
Describes a single nucleotide.
"""
abstract type SingleNucleotide <: AbstractDescriptor end
"""
Describes a site inbetween two nucleotides, e.g. a cleavage site ("1^2"). `position` points to the nucleotide to the left of the site ("1" in the previous example).
"""
abstract type BetweenNucleotides <: AbstractDescriptor end
"""
Describes a normal locus, e.g. "1..1000"
"""
abstract type ClosedSpan <: AbstractDescriptor end
abstract type OpenRightSpan <: AbstractDescriptor end
abstract type OpenLeftSpan <: AbstractDescriptor end
abstract type OpenSpan <: AbstractDescriptor end
"""
Not yet implemented.
Describes a remote entry, e.g. exons found on another chromosome.
"""
abstract type Remote end

struct PointLocus{T} <: AbstractLocus where {T <: AbstractDescriptor}
    position::Int
    descriptor::Type{T}
end

struct SpanLocus{T} <: AbstractLocus where {T <: AbstractDescriptor}
    position::UnitRange{Int}
    descriptor::Type{T}
end

struct Join{L <: AbstractLocus} <: AbstractLocus
    loc::Vector{L}
end

struct Order{L <: AbstractLocus} <: AbstractLocus
    loc::Vector{L}
end

struct Complement{L <: Union{PointLocus{T} where T, SpanLocus{T} where T, Join, Order}} <: AbstractLocus
    loc::L
end

Locus() = PointLocus(1, SingleNucleotide)

ClosedSpan(p) = SpanLocus(p, ClosedSpan)
OpenRightSpan(p) = SpanLocus(p, OpenRightSpan)
OpenLeftSpan(p) = SpanLocus(p, OpenLeftSpan)
OpenSpan(p) = SpanLocus(p, OpenSpan)

SingleNucleotide(p) = PointLocus(p, SingleNucleotide)
BetweenNucleotides(p) = PointLocus(p, BetweenNucleotides)

Complement(loc::Complement) = loc.loc

Base.convert(::Type{ClosedSpan}, p::UnitRange{Int}) = ClosedSpan(p)
Base.convert(::Type{Complement{L}}, p::UnitRange{Int}) where L <: AbstractLocus = Complement(L(p))
function Base.convert(::Type{ClosedSpan}, p::StepRange{Int, Int})
    if p.step == -1
        return Complement(ClosedSpan(p.stop:p.start))
    elseif p.step == 1
        return ClosedSpan(p.start:p.stop)
    end
    throw(DomainError(1, "`p` must have a step of 1 or -1"))
end


Base.show(io::IO, locus::PointLocus{SingleNucleotide}) = print(io, string(locus.position))
Base.show(io::IO, locus::PointLocus{BetweenNucleotides}) = print(io, string(locus.position, "^", locus.position + 1))
Base.show(io::IO, locus::SpanLocus{ClosedSpan}) = print(io, string(locus.position.start, "..", locus.position.stop))
Base.show(io::IO, locus::SpanLocus{OpenSpan}) = print(io, string("<", locus.position.start, "..", ">", locus.position.stop))
Base.show(io::IO, locus::SpanLocus{OpenLeftSpan}) = print(io, string("<", locus.position.start, "..", locus.position.stop))
Base.show(io::IO, locus::SpanLocus{OpenRightSpan}) = print(io, string(locus.position.start, "..", ">", locus.position.stop))
Base.show(io::IO, locus::Join) = print(io, string("join(", join(locus.loc, ","), ")"))
Base.show(io::IO, locus::Order) = print(io, string("order(", join(locus.loc, ","), ")"))
function Base.show(io::IO, locus::Complement)
    print(io, "complement(", string(locus.loc), ")")
end


Base.:(==)(loc1::SpanLocus{ClosedSpan}, loc2::SpanLocus{ClosedSpan}) = loc1.position == loc2.position
Base.:(==)(loc1::PointLocus{SingleNucleotide}, loc2::PointLocus{SingleNucleotide}) = loc1.position == loc2.position
Base.:(==)(loc1::PointLocus{BetweenNucleotides}, loc2::PointLocus{BetweenNucleotides}) = loc1.position == loc2.position
Base.:(==)(loc1::Join{T}, loc2::Join{T}) where {T <: AbstractLocus} = (length(loc1.loc) == length(loc2.loc)) && all(pair -> pair[1] == pair[2], zip(loc1.loc, loc2.loc))
Base.:(==)(loc1::Order{T}, loc2::Order{T}) where {T <: AbstractLocus} = (length(loc1.loc) == length(loc2.loc)) && all(pair -> pair[1] == pair[2], zip(loc1.loc, loc2.loc))
Base.:(==)(loc1::Complement{T}, loc2::Complement{T}) where {T <: AbstractLocus} = loc1.loc == loc2.loc
Base.:(==)(loc1::AbstractLocus, loc2::AbstractLocus) = false

Base.iterate(loc::Union{PointLocus, SpanLocus}) = (loc, nothing)
Base.iterate(loc::Union{PointLocus, SpanLocus}, ::Any) = nothing
Base.iterate(loc::Union{Join, Order}, i = 1) = i > length(loc.loc) ? nothing : (loc.loc[i], i+1)
Base.iterate(loc::Complement{T}) where {T <: Union{PointLocus, SpanLocus}} = (loc, 1)
Base.iterate(loc::Complement{T}, ::Any) where {T <: Union{PointLocus, SpanLocus}} = nothing
Base.iterate(loc::Complement{T}, i = 1) where {T <: Union{Join, Order}} = i > length(loc.loc.loc) ? nothing : (Complement(loc.loc.loc[end-i+1]), i+1)

Base.IteratorSize(loc::AbstractLocus) = Base.SizeUnknown()

function Base.isless(l1::AbstractLocus, l2::AbstractLocus)
    if (l1.start < l2.start) || (l1.start == l2.start && l1.stop > l2.stop)
        return true
    end
    return false
end

"""
    eachposition(loc::AbstractLocus)

Return an object that iterates over each position in the locus in the specified order. Returns `nothing` for `PointLocus{BetweenNucleotides}`.

```julia
julia> eachposition(Locus("join(1..3,complement(7..9))"))
[1,2,3,9,8,7]
"""
eachposition(loc::PointLocus{BetweenNucleotides}) = nothing
eachposition(loc::PointLocus{SingleNucleotide}) = loc.position
eachposition(loc::SpanLocus) = loc.position
eachposition(loc::Union{Join, Order}) = Iterators.flatten(map(eachposition, loc.loc))
eachposition(loc::Complement{T}) where {T <: Union{PointLocus, SpanLocus}} = Iterators.reverse(loc.position)
eachposition(loc::Complement) = Iterators.reverse(eachposition(loc.loc))


"""
    Locus(s::AbstractString)

Parse `s` as a GenBank location descriptor and return an appropriate `AbstractLocus`.
"""
function Locus(s::T) where T <: AbstractString
    ### Complement
    if occursin(r"^complement", s)
        return Complement(Locus(s[12:end-1]))
    end
    ### Join
    if occursin(r"^join", s)
        S = split(s[6:end-1], ',')
        return Join([Locus(z) for z in S])
    end
    if occursin(r"^order", s)
        S = split(s[7:end-1], ',')
        return Order([Locus(z) for z in S])
    end
    ### SingleNucleotide
    occursin(r"^\d+$", s) && return PointLocus(parse(Int, s), SingleNucleotide)
    ### BetweenNucleotides
    m = match(r"^(\d+)\^", s)
    if !isnothing(m)
        return PointLocus(parse(Int, m[1]), BetweenNucleotides)
    end
    ### Span
    m = match(r"^(<?)(\d+)\.\.(>?)(\d+)$", s)
    if !isnothing(m)
        p = parse(Int, m[2]):parse(Int, m[4])
        type = if all(isempty, m.captures[[1,3]])
            ClosedSpan
        elseif all(!isempty, m.captures[[1,3]])
            OpenSpan
        elseif !isempty(m.captures[1])
            OpenLeftSpan
        else
            OpenRightSpan
        end
        return SpanLocus(p, type)
    end
    m = match(r"^(<?)(>?)(\d+)$", s)
    if !isnothing(m)
        p = parse(Int, m[3])
        type = if !isempty(m[1]) && !isempty(m[3])
            OpenSpan
        elseif !isempty(m[1])
            OpenLeftSpan
        else
            OpenRightSpan
        end
        return SpanLocus(p:p, type)
    end
    @warn "Couldn't parse location descriptor \"$s\""
    return nothing
end


"""
    shift(loc, p)

Shift the position of `Locus` loc by `p`.
"""
function shift(loc, p) 
    
end