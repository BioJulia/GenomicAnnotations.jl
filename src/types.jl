"""
Abstract type used to represent genes on a `Chromosome`. Subtypes must implement:
 * parent(gene)
 * index(gene)
 * locus(gene)
 * addgene!(chr, feature, locus)

This mainly exists to allow `Chromosome`s and `Gene`s to refer to eachother,
and may be removed in the future.
"""
abstract type AbstractGene end


"""
Struct for storing information on genomic locations. `strand` can be '+', '-',
or '.' when the strand is irrelevant. `order` is used to store discontiguous
sequences, indicated in the GenBank file with the order() and join() operators.
"""
struct Locus
    position::UnitRange{Int}
    strand::Char
    complete_left::Bool
    complete_right::Bool
    order::Vector{UnitRange{Int}}
    join::Bool
end


"""
    Locus()
    Locus(position::UnitRange{Int})
    Locus(position::UnitRange{Int}, strand::Char)

"""
Locus() = Locus(1:1, '.', true, true, UnitRange{Int}[], false)
Locus(position::UnitRange{Int}) = Locus(position, '.', true, true, UnitRange{Int}[], false)
Locus(position::UnitRange{Int}, strand::Char) = Locus(position, strand, true, true, UnitRange{Int}[], false)
Locus(position::UnitRange{Int}, strand::Char, complete_left, complete_right) = Locus(position, strand, complete_left, complete_right, UnitRange{Int}[], false)

Base.convert(::Type{Locus}, x::UnitRange{Int}) = Locus(x)
function Base.convert(::Type{Locus}, x::StepRange{Int, Int})
    if x.step == -1
        return Locus(x.stop:x.start, '-')
    elseif x.step == 1
        return Locus(x.start:x.stop)
    end
    throw(DomainError(x, "`x` must have a step of 1 or -1"))
end


"""
Struct for storing annotations for a single chromosome, plasmid, contig, etc.
Contains five fields: `name`, `sequence`, `header`, `genes`, and `genedata`.
Annotations are stored as a `DataFrame` in `genedata`, but can be accessed
more easily through `genes` using the API provided in this module.
"""
mutable struct Chromosome{G <: AbstractGene}
    name::String
    sequence::LongDNASeq
    header::String
    genes::Vector{G}
    genedata::DataFrame
    circular::Bool
    Chromosome{G}(name, sequence, header, genes::Vector{G}, genedata, circular) where G = new(name, sequence, header, genes, genedata, circular)
end


Chromosome(args...) = Chromosome{Gene}(args...)


struct Gene <: AbstractGene
    parent::Chromosome{Gene}
    index::UInt
    locus::Locus
    feature::Symbol
end


Chromosome{Gene}() = Chromosome{Gene}("", dna"", "", Gene[], DataFrame(), false)


iscircular(c::Chromosome{T}) where {T<:AbstractGene} = c.circular


"""
    addgene!(chr::Chromosome, feature, locus; kw...)

Add gene to `chr`. `locus` can be a `Locus`, a UnitRange, or a StepRange (for
decreasing ranges, which will be annotated on the complementary strand).

# Example
```julia
addgene!(chr, "CDS", 1:756;
    locus_tag = "gene0001",
    product = "Chromosomal replication initiator protein dnaA")
```
"""
function addgene!(chr::Chromosome{Gene}, feature, locus; kw...)
    locus = convert(Locus, locus)
    push!(chr.genedata, fill(missing, size(chr.genedata, 2)))
    index = UInt32(length(chr.genes) + 1)
    gene = Gene(chr, index, locus, feature)
    for (k, v) in kw
        Base.setproperty!(gene, k, v)
    end
    push!(chr.genes, gene)
    return gene
end


"""
    delete!(gene::AbstractGene)

Delete `gene` from `parent(gene)`. Warning: does not work when broadcasted! Use
`delete!(::AbstractVector{Gene}) instead.`
"""
function Base.delete!(gene::AbstractGene)
    i = index(gene)
    deleterows!(parent(gene).genedata, i)
    deleteat!(parent(gene).genes, lastindex(parent(gene).genes))
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
    indices = index(gene)
    DataFrames.deleterows!(parent(genes).genedata, Int.(indices))
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
function pushproperty!(gene::AbstractGene, qualifier::Symbol, value::T; forceany = true) where T
    if hasproperty(parent(gene).genedata, qualifier)
        C = eltype(parent(gene).genedata[!, qualifier])
        if T <: C
            if ismissing(parent(gene).genedata[index(gene), qualifier])
                parent(gene).genedata[index(gene), qualifier] = value
            else
                parent(gene).genedata[!, qualifier] = vectorise(parent(gene).genedata[!, qualifier])
                push!(parent(gene).genedata[index(gene), qualifier], value)
            end
        elseif Vector{T} <: C
            if ismissing(parent(gene).genedata[index(gene), qualifier])
                parent(gene).genedata[index(gene), qualifier] = [value]
            else
                push!(parent(gene).genedata[index(gene), qualifier], value)
            end
        elseif forceany && !(C <: AbstractVector)
            if ismissing(parent(gene).genedata[index(gene), qualifier])
                parent(gene).genedata[!, qualifier] = convert(Vector{Any}, parent(gene).genedata[!, qualifier])
                parent(gene).genedata[index(gene), qualifier] = value
            else
                parent(gene).genedata[!, qualifier] = vectorise(convert(Vector{Any}, parent(gene).genedata[! ,qualifier]))
                push!(parent(gene).genedata[index(gene), qualifier], value)
            end
        elseif forceany && C <: AbstractVector
            if ismissing(parent(gene).genedata[index(gene), qualifier])
                parent(gene).genedata[!, qualifier] = convert(Vector{Any}, parent(gene).genedata[!, qualifier])
                parent(gene).genedata[index(gene), qualifier] = [value]
            else
                @error "This shouldn't happen"
            end
        else
            @error "Tried to add a '$T' to '$qualifier::$(typeof(parent(gene).genedata[!, qualifier]))'"
        end
    else
        s = size(parent(gene).genedata, 1)
        parent(gene).genedata[!, qualifier] = Vector{Union{Missing, T}}(missing, s)
        parent(gene).genedata[index(gene), qualifier] = value
    end
    return value
end


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
    sequence(gene::AbstractGene; translate = false)

Return genomic sequence for `gene`. If `translate` is `true`, the sequence will be translated to a `LongAminoAcidSeq`, excluding the stop, otherwise it will be returned as a `LongDNASeq` (including the stop codon).
```
"""
function sequence(gene::AbstractGene; translate = false)
    if locus(gene).strand == '-'
        s = reverse_complement(parent(gene).sequence[locus(gene).position])
    else
        s = parent(gene).sequence[locus(gene).position]
    end
    translate ? BioSequences.translate(s)[1:end-1] : s
end


Base.length(gene::AbstractGene) = length(locus(gene))
Base.length(locus::Locus) = length(locus.position)


"""
    iscomplement(gene)

Return `true` if `locus(gene).compliment == '-'`, otherwise return `false`.
"""
iscomplement(gene::AbstractGene) = locus(gene).strand == '-'


"""
    iscomplete(gene)

Return `true` if `gene` is a complete gene, i.e. not a pseudo gene or partial.
"""
iscomplete(gene::AbstractGene) = !any(get(gene, :pseudo, false)) && !any(get(gene, :ribosomal_slippage, false)) && locus(gene).complete_right && locus(gene).complete_left


function appendstring(field, v::Bool)
    return "\n" * join(fill(' ', 21)) * "/$field"
end
function appendstring(field, v::Union{Number, Symbol})
    return "\n" * join(fill(' ', 21)) * "/$field=" * string(v)
end
function appendstring(field, v)
	v = _multiline(v, field)
    v = replace(v, "\n" => "\n" * join(fill(' ', 21)))
    return "\n" * join(fill(' ', 21)) * "/$field=\"" * v * "\""
end


function _findbreak(v, t)
	i = t == 59 ?
		findlast('\n', v) :
		1
	j = i+t-1
	# j >= lastindex(v) && return nothing
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
	_multiline(v, s::Symbol)

Return a `String` with newlines and spaces added so that it conforms to the GenBank line width.
"""
function _multiline(v, s)
	v = string(v)
	s = string(s)
	if length(v) + length(s) > 55 && !occursin('\n', v)
		b = _findbreak(v, 55 - length(s))
		v = v[1:b.start] * "\n" * v[b.stop:end]
		b = _findbreak(v, 59)
		while !isnothing(b) && b.start < length(v)
			v = v[1:b.start] * "\n" * v[b.stop:end]
			b = _findbreak(v, 59)
		end
	end
	return v
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


function Base.show(io::IO, chr::Chromosome)
    s = "Chromosome '$(chr.name)' ($(length(chr.sequence)) bp) with $(length(chr.genes)) annotations\n"
    print(io, s)
end


function Base.show(io::IO, genes::Vector{AbstractGene})
    print(join(genes, "\n"))
end


function Base.show(io::IO, locus::Locus)
    s = ""
    locus.strand == '-'     && (s *= "complement(")
    !locus.complete_left    && (s *= ">")
    if length(locus.order) > 0
        if locus.join
            s *= "join(" * join([join((r.start, r.stop), "..") for r in locus.order], ",") * ")"
        else
            s *= "order(" * join([join((r.start, r.stop), "..") for r in locus.order], ",") * ")"
        end
    else
        s *= join((locus.position.start, locus.position.stop), "..")
    end
    !locus.complete_right   && (s *= "<")
    locus.strand == '-'     && (s *= ")")
    print(io, s)
end


Base.copy(chr::Chromosome) = deepcopy(chr)
Base.copy(gene::AbstractGene) = deepcopy(gene)


function Base.isless(l1::Locus, l2::Locus)
    if (l1.position.start < l2.position.start) || (l1.position.start == l2.position.start && l1.position.stop > l2.position.stop)
        return true
    end
    return false
end
Base.isless(g1::AbstractGene, g2::AbstractGene) = ((locus(g1) == locus(g2)) && (feature(g1) == "gene" && feature(g2) != "gene")) || (locus(g1) < locus(g2))


function Base.:(==)(x::Locus, y::Locus)
    (x.position == y.position) && (x.strand == y.strand) && (x.complete_left == y.complete_left) && (x.complete_right == y.complete_right) && (x.order == y.order)
end

Base.in(loc::Locus, r::UnitRange{Int}) = loc.position.start in r && loc.position.stop in r
Base.in(loc1::Locus, loc2::Locus) = loc1 in loc2.position
Base.intersect(loc1::Locus, loc2::Locus) = intersect(loc1.position, loc2.position)
Base.iterate(loc::Locus) = iterate(loc.position)

index(g::Gene) = getfield(g, :index)
locus(g::Gene) = getfield(g, :locus)
feature(g::Gene) = getfield(g, :feature)
Base.parent(g::Gene) = getfield(g, :parent)
function Base.parent(gs::AbstractVector{G}) where {G <: AbstractGene}
    @assert !isempty(gs) "Trying to get parent of empty Array{Gene}"
    getfield(gs[1], :parent)
end
