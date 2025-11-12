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


mutable struct Gene <: AbstractGene
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
            gd[!, qualifier] = Union{Missing, T}[value] :
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


"""
    sequence(record::Record)

Return the nucleotide sequence of `record`.
"""
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
    iscomplete(gene::AbstractGene)
    iscomplete(loc::AbstractLocus)

Return `true` if `gene` is a complete gene, i.e. not a pseudo gene or partial. Compound loci (i.e. those with join/order) are considered incomplete if any element is incomplete.
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

"""
    iscompound(gene::AbstractGene)
    iscompound(loc::AbstractLocus)

Return `true` if `loc` is a `Join`, `Order`, or either wrapping in `Complement`.
"""
iscompound(gene::AbstractGene) = iscompound(locus(gene))
iscompound(locus::Union{SpanLocus, PointLocus}) = false
iscompound(locus::Union{Join, Order}) = true
iscompound(locus::Complement{T}) where T <: AbstractLocus = iscompound(locus.loc)

function appendstring(field, v::Bool)
    return "\n" * join(fill(' ', 21)) * "/$field"
end
function appendstring(field, v::Union{Number, Symbol})
    return "\n" * join(fill(' ', 21)) * "/$field=" * string(v)
end
function appendstring(field, v)
    v = if v[[1,end]] == ['(', ')']
        string(v)
    else
        string("\"", v, "\"")
    end
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

Base.show(io::IO, mime::MIME"text/plain", genes::Array{G, 1}) where {G <: AbstractGene} = Base.show(io, "Vector{$G} containing $(length(genes)) annotations")
function Base.show(io::IO, genes::Array{G, 1}) where {G <: AbstractGene}
    for gene in genes
        show(io, gene)
    end
end

function Base.show(io::IO, chr::Record)
    s = "Chromosome '$(chr.name)' ($(length(chr.sequence)) bp) with $(length(chr.genes)) annotations\n"
    print(io, s)
end

Base.copy(chr::Record) = deepcopy(chr)
Base.copy(gene::AbstractGene) = deepcopy(gene)


"""
    index(g::Gene)

Return the index of `g`, i.e. its row number in `parent(g).genedata`.
"""
index(g::Gene) = getfield(g, :index)

"""
    locus(g::Gene)

Return the `AbstractLocus` of `g`.
"""
locus(g::Gene) = getfield(g, :locus)
Base.position(g::Gene) = getfield(g, :locus).position

"""
    feature(g::Gene)

Return the feature type (i.e. gene, CDS, tRNA, etc.) of `g`.
"""
feature(g::Gene) = getfield(g, :feature)

"""
    parent(g::Gene)
    parent(gs::AbstractVector{Gene})

Return the parent `Record` of `g`. Errors for `AbstractVector{Gene}`s if the genes do not come from the same parent.
"""
Base.parent(g::Gene) = getfield(g, :parent)
function Base.parent(gs::AbstractVector{G}) where {G <: AbstractGene}
    @assert length(unique(map(parent, gs))) == 1 "Trying to get parent of a Vector{Gene} with genes from multiple `Record`s."
    @assert !isempty(gs) "Trying to get parent of empty Array{Gene}"
    getfield(gs[1], :parent)
end

"""
    genedata(g::Gene)

Return the `DataFrameRow` where the data for `g` is stored. See [`attributes`](@ref) for a slower but potentially more convenient alternative.
"""
genedata(g::Gene) = parent(g).genedata[index(g), :]

"""
    attributes(g::Gene)

Return an immutable NamedTuple containing copies of all annotated attributes of `g`. Missing attributes are excluded. See [`genedata`](@ref) for a non-allocating way to access the gene data directly.
"""
attributes(g::Gene) = (; filter(kw -> !ismissing(kw[2]), collect(pairs(genedata(g))))...)

"""
    feature!(g::Gene, f::Symbol)

Change the feature of `g` to `f`, returning a new instance of `Gene`. Since `Gene`s are immutable, `feature!` only mutates the parent of `g` and not `g` itself. Thus, in the first example below the original unmodified `g` is printed, not the updated version:

```julia
# This will not work as expected:
for source in @genes(chr, source)
    feature!(source, :region)
    println(source)
end

# But this will:
for source in @genes(chr, source)
    source = feature!(source, :region)
    println(source)
```
"""
function feature!(g::Gene, f::Symbol)
    if f == feature(g)
        return g
    else
        parent(g).genes[index(g)] = Gene(parent(g), index(g), locus(g), f)
    end
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

Base.isless(g1::AbstractGene, g2::AbstractGene) = ((locus(g1) == locus(g2)) && (feature(g1) == "gene" && feature(g2) != "gene")) || (locus(g1) < locus(g2))

function shift!(gene, p)
    setfield!(gene, :locus, shift(locus(gene), p))
    return gene
end