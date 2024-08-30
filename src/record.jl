"""
Abstract type used to represent genes on a `Record`. Subtypes must implement:
 * parent(gene)
 * index(gene)
 * locus(gene)
 * addgene!(chr, feature, locus)

This mainly exists to allow `Record`s and `Gene`s to refer to eachother,
and may be removed in the future.
"""
abstract type AbstractGene end


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


"""
Struct for storing annotations for a single chromosome, plasmid, contig, etc.
Contains five fields: `name`, `sequence`, `header`, `genes`, and `genedata`.
Annotations are stored as a `DataFrame` in `genedata`, but can be accessed
more easily through `genes` using the API provided in this module.
"""
mutable struct Record{G <: AbstractGene}
    name::String
    sequence::LongDNA{4}
    header::String
    genes::Vector{G}
    genedata::DataFrame
    circular::Bool
    Record{G}(name, sequence, header, genes::Vector{G}, genedata, circular) where G = new(name, sequence, header, genes, genedata, circular)
end


Record(args...) = Record{Gene}(args...)


struct Gene <: AbstractGene
    parent::Record{Gene}
    index::UInt
    locus::AbstractLocus
    feature::Symbol
end


Record{Gene}() = Record{Gene}("", dna"", "", Gene[], DataFrame(), false)


iscircular(c::Record{T}) where {T<:AbstractGene} = c.circular


"""
    addgene!(chr::Record, feature, locus; kw...)

Add gene to `chr`. `locus` can be an `AbstractLocus`, a String, a UnitRange, or a StepRange (for
decreasing ranges, which will be annotated on the complementary strand).

# Example
```julia
addgene!(chr, "CDS", 1:756;
    locus_tag = "gene0001",
    product = "Chromosomal replication initiator protein dnaA")
```
"""
function addgene!(chr::Record{Gene}, feature, loc::AbstractLocus, geneindex; kw...)
    push!(chr.genedata, fill(missing, size(chr.genedata, 2)))
    gene = Gene(chr, geneindex, loc, feature)
    for (k, v) in kw
        pushproperty!(chr, geneindex, k, v)
    end
    push!(chr.genes, gene)
    return gene
end
addgene!(chr::Record{Gene}, feature, loc::AbstractString, geneindex; kw...) = addgene!(chr, feature, Locus(loc), geneindex; kw...)
addgene!(chr::Record{Gene}, feature, loc; kw...) = addgene!(chr, feature, loc, UInt32(length(chr.genes) + 1); kw...)


"""
    delete!(gene::AbstractGene)

Delete `gene` from `parent(gene)`. Warning: does not work when broadcasted! Use
`delete!(::AbstractVector{Gene}) instead.`
"""
function Base.delete!(gene::AbstractGene)
    i = index(gene)
    deleteat!(parent(gene).genes, i)
    nothing
end


"""
    delete!(genes::AbstractArray{Gene, 1})

Delete all genes in `genes` from `parent(genes[1])`.

# Example
```julia
delete!(@genes(chr, length(gene) <= 60))
```
"""
function Base.delete!(genes::AbstractArray{Gene, 1})
    indices = index.(genes)
    DataFrames.delete!(parent(genes).genedata, Int.(indices))
    lastindices = length(parent(genes).genes) - length(indices) + 1 : length(parent(genes).genes)
    deleteat!(parent(genes).genes, lastindices)
    nothing
end

function Base.delete!(genes::Vector{Gene})
    indices = index.(genes)
    DataFrames.delete!(parent(genes).genedata, Int.(indices))
    lastindices = length(parent(genes).genes) - length(indices) + 1 : length(parent(genes).genes)
    deleteat!(parent(genes).genes, lastindices)
    nothing
end

function Base.propertynames(gene::G) where {G <: AbstractGene}
    propertynames(parent(gene).genedata)
end


"""
    vectorise(A::AbstractArray{Union{Missing, T}, 1}) where T

Convert an array of type Vector{Union{Missing, T}} to Vector{Union{Missing, Vector{T}}}.
"""
function vectorise(A::AbstractArray{Union{Missing, T}, 1}) where T
    B = Union{Missing, Vector{T}}[]
    for a in A
        if ismissing(a)
            push!(B, missing)
        else
            push!(B, [a])
        end
    end
    return B
end


"""
    pushproperty!(gene::AbstractGene, qualifier::Symbol, value::T)

Add a property to `gene`, similarly to `Base.setproperty!(::gene)`, but if the
property is not missing in `gene`, it will be transformed to store a vector
instead of overwriting existing data.

```julia
julia> eltype(chr.genedata[!, :EC_number])
Union{Missing,String}

julia> chr.genes[1].EC_number = "EC:1.2.3.4"
"EC:1.2.3.4"

julia> pushproperty!(chr.genes[1], :EC_number, "EC:4.3.2.1"); chr.genes[1].EC_number
2-element Array{String,1}:
 "EC:1.2.3.4"
 "EC:4.3.2.1"

julia> eltype(chr.genedata[!, :EC_number])
Union{Missing, Array{String,1}}
```
"""
function pushproperty!(chr::Record, ind::Integer, qualifier::Symbol, value::T; forceany = true) where T
    gd = chr.genedata
    if hasproperty(gd, qualifier)
        C = eltype(gd[!, qualifier])
        if T <: C
            if ismissing(gd[ind, qualifier])
                gd[ind, qualifier] = value
            else
                gd[!, qualifier] = vectorise(gd[!, qualifier])
                push!(gd[ind, qualifier], value)
            end
        elseif Vector{T} <: C
            if ismissing(gd[ind, qualifier])
                gd[ind, qualifier] = [value]
            else
                push!(gd[ind, qualifier], value)
            end
        elseif forceany && !(C <: AbstractVector)
            if ismissing(gd[ind, qualifier])
                gd[!, qualifier] = convert(Vector{Any}, gd[!, qualifier])
                gd[ind, qualifier] = value
            else
                gd[!, qualifier] = vectorise(convert(Vector{Any}, gd[!, qualifier]))
                push!(gd[ind, qualifier], value)
            end
        elseif forceany && C <: AbstractVector
            if ismissing(gd[ind, qualifier])
                gd[!, qualifier] = convert(Vector{Any}, gd[!, qualifier])
                gd[ind, qualifier] = [value]
            else
                @error "This shouldn't happen"
            end
        else
            @error "Tried to add a '$T' to '$qualifier::$(typeof(gd[!, qualifier]))'"
        end
    else
        s = size(gd, 1)
        isempty(names(gd)) ?
            gd[!, qualifier] = Union{Missing, typeof(value)}[value] :
            gd[!, qualifier] = missings(T, s)
        gd[ind, qualifier] = value
    end
    return value
end
pushproperty!(gene::AbstractGene, qualifier::Symbol, value::T; forceany = true) where T = pushproperty!(parent(gene), index(gene), qualifier, value; forceany = forceany)

"""
    get(g::AbstractGene, key, default)

Retrieve `key` from `g`. If `key` is missing, return `default`.
"""
function Base.get(g::AbstractGene, key, default)
    p = getproperty(g, key)
    if ismissing(p)
        return default
    else
        return p
    end
end


"""
    sequence(gene::AbstractGene; translate = false, preserve_alternate_start = false)

Return genomic sequence for `gene`. If `translate` is `true`, the sequence will be translated to a `LongAA`, excluding the stop, otherwise it will be returned as a `LongDNA{4}` (including the stop codon). If `preserve_alternate_start` is set to false, alternate start codons will be assumed to code for methionine.
```
"""
sequence(gene::AbstractGene; kw...) = sequence(parent(gene).sequence, locus(gene); kw...)
function sequence(chrseq, loc::AbstractLocus; translate = false, preserve_alternate_start = false)
    seq = _sequence(chrseq, loc)
    if translate
        if preserve_alternate_start
            return BioSequences.translate(seq[1:end-3])
        end
        return aa"M" * BioSequences.translate(seq[4:end-3])
    end
    return seq
end

_sequence(chrseq, loc::SpanLocus) = @view(chrseq[loc.position])
_sequence(chrseq, loc::PointLocus{SingleNucleotide}) = @view(chrseq[loc.position:loc.position])
_sequence(chrseq, loc::PointLocus{BetweenNucleotides}) = @view(chrseq[loc.position:(loc.position + 1)])
_sequence(chrseq, loc::Complement) = reverse_complement(_sequence(chrseq, loc.loc))
_sequence(chrseq, loci::Complement{T}) where T <: Order = map(reverse_complement, reverse(_sequence(chrseq, loci.loc)))
_sequence(chrseq::T, loci::Join) where T = *(map(seq -> convert(T, seq), [_sequence(chrseq, loc) for loc in loci.loc])...)
_sequence(chrseq, loci::Order) = [_sequence(chrseq, loc) for loc in loci.loc]


Base.length(gene::AbstractGene) = length(locus(gene))
Base.length(locus::PointLocus) = 1
Base.length(locus::SpanLocus) = length(locus.position)
Base.length(locus::Union{Join, Order}) = sum(map(length, locus.loc))
Base.length(locus::Complement) = length(locus.loc)


function sequence(record::Record)
    record.sequence
end


"""
    iscomplement(gene::AbstractGene)
    iscomplement(loc::AbstractLocus)

Return `true` if `loc`/`locus(gene)` is a `Compliment`, or if it is a `Join` or `Order` whose elements are all `Complement`s. Otherwise return `false`.
"""
iscomplement(gene::AbstractGene) = iscomplement(locus(gene))
iscomplement(loc::Union{Join, Order}) = all(iscomplement, loc.loc)
iscomplement(loc::Complement) = true
iscomplement(loc::AbstractLocus) = false


"""
    iscomplete(gene)

Return `true` if `gene` is a complete gene, i.e. not a pseudo gene or partial.
"""
iscomplete(gene::AbstractGene) = iscomplete(locus(gene)) && !any(get(gene, :pseudo, false)) && !any(get(gene, :partial, false))
iscomplete(locus::SpanLocus{ClosedSpan}) = true
iscomplete(locus::SpanLocus{OpenSpan}) = false
iscomplete(locus::SpanLocus{OpenRightSpan}) = false
iscomplete(locus::SpanLocus{OpenLeftSpan}) = false
iscomplete(locus::PointLocus{T}) where T = true
iscomplete(locus::Complement) = iscomplete(locus.loc)
iscomplete(loci::Join) = all(iscomplete, loci.loc)
iscomplete(loci::Order) = all(iscomplete, loci.loc)


function appendstring(field, v::Bool)
    return "\n" * join(fill(' ', 21)) * "/$field"
end
function appendstring(field, v::Union{Number, Symbol})
    return "\n" * join(fill(' ', 21)) * "/$field=" * string(v)
end
function appendstring(field, v)
    v = string("\"", v, "\"")
    v = multiline(v, field)
    v = replace(v, "\n" => "\n" * join(fill(' ', 21)))
    return "\n" * join(fill(' ', 21)) * "/$field=" * v
end


function Base.show(io::IO, gene::AbstractGene)
    buf = IOBuffer()
    print(buf, "     " * rpad(string(feature(gene)), 16, ' '))
    print(buf, locus(gene))
    for field in names(parent(gene).genedata)
        field in [:feature, :locus] && continue
        v = parent(gene).genedata[index(gene), field]
        if !ismissing(v)
            if v isa AbstractVector
                for i in eachindex(v)
                    print(buf, appendstring(field, v[i]))
                end
            else
                print(buf, appendstring(field, v))
            end
        end
    end
    println(io, String(take!(buf)))
end


Base.show(io::IO, mime::MIME"text/plain", genes::Array{G, 1}) where {G <: AbstractGene} = Base.show(io, genes)
function Base.show(io::IO, genes::Array{G, 1}) where {G <: AbstractGene}
    for gene in genes
        show(io, gene)
    end
end


function Base.show(io::IO, chr::Record)
    s = "Chromosome '$(chr.name)' ($(length(chr.sequence)) bp) with $(length(chr.genes)) annotations\n"
    print(io, s)
end


function Base.show(io::IO, genes::Vector{AbstractGene})
    print(join(genes, "\n"))
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

Base.copy(chr::Record) = deepcopy(chr)
Base.copy(gene::AbstractGene) = deepcopy(gene)


function Base.isless(l1::AbstractLocus, l2::AbstractLocus)
    if (l1.start < l2.start) || (l1.start == l2.start && l1.stop > l2.stop)
        return true
    end
    return false
end
Base.isless(g1::AbstractGene, g2::AbstractGene) = ((locus(g1) == locus(g2)) && (feature(g1) == "gene" && feature(g2) != "gene")) || (locus(g1) < locus(g2))


Base.:(==)(loc1::SpanLocus{ClosedSpan}, loc2::SpanLocus{ClosedSpan}) = loc1.position == loc2.position
Base.:(==)(loc1::PointLocus{SingleNucleotide}, loc2::PointLocus{SingleNucleotide}) = loc1.position == loc2.position
Base.:(==)(loc1::Join{T}, loc2::Join{T}) where {T <: AbstractLocus} = (length(loc1.loc) == length(loc2.loc)) && all(pair -> pair[1] == pair[2], zip(loc1.loc, loc2.loc))
Base.:(==)(loc1::Order{T}, loc2::Order{T}) where {T <: AbstractLocus} = (length(loc1.loc) == length(loc2.loc)) && all(pair -> pair[1] == pair[2], zip(loc1.loc, loc2.loc))
Base.:(==)(loc1::Complement{T}, loc2::Complement{T}) where {T <: AbstractLocus} = loc1.loc == loc2.loc
Base.:(==)(loc1::AbstractLocus, loc2::AbstractLocus) = false

Base.in(loc::PointLocus{SingleNucleotide}, r::UnitRange) = loc.position in r
Base.in(loc::PointLocus{BetweenNucleotides}, r::UnitRange) = loc.position in r[1:end-1]
Base.in(loc::SpanLocus, r::UnitRange) = loc.position in r
Base.in(loc::Join, r::UnitRange) = all(in.(loc.loc, r))
Base.in(loc::Order, r::UnitRange) = all(in.(loc.loc, r))
Base.in(loc::Complement, r::UnitRange) = all(in.(loc.loc, r))

# Base.intersect(loc1::PointLocus{SingleNucleotide}, loc2::A)
# Base.intersect(loc1::Locus, loc2::Locus) = intersect(loc1.position, loc2.position)

Base.iterate(loc::PointLocus{SingleNucleotide}) = iterate(loc.position)
Base.iterate(loc::SpanLocus{T}) where T = iterate(loc.position)
Base.iterate(loc::AbstractLocus) = iterate(union(x.loc for x in loc.loc))

index(g::Gene) = getfield(g, :index)
locus(g::Gene) = getfield(g, :locus)
Base.position(g::Gene) = getfield(g, :locus).position
feature(g::Gene) = getfield(g, :feature)
Base.parent(g::Gene) = getfield(g, :parent)
function Base.parent(gs::AbstractVector{G}) where {G <: AbstractGene}
    @assert !isempty(gs) "Trying to get parent of empty Array{Gene}"
    getfield(gs[1], :parent)
end


"""
    locus!(gene::AbstractGene, loc)
    locus!(gene::AbstractGene, loc::AbstractLocus)

Replace `gene` with a new `Gene` with `loc` as its `Locus`. If `loc` is not an `AbstractLocus`, it is parsed with `Locus(loc)`.
"""
locus!(gene::AbstractGene, x) = locus!(gene, Locus(x))
function locus!(gene::AbstractGene, loc::AbstractLocus)
    chr = parent(gene)
    newgene = Gene(chr, index(gene), loc, feature(gene))
    setindex!(chr.genes, newgene, index(gene))
    newgene
end


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
    m = match(r"^(<?)(\d+)..(>?)(\d+)$", s)
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