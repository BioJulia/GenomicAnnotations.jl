abstract type AbstractGene end


"""
Struct for storing information on genomic locations. `strand` can be '+', '-',
or '.' when the strand is irrelevant.
"""
struct Locus
    position::UnitRange{Int}
    strand::Char
    complete_left::Bool
    complete_right::Bool
    order::Vector{UnitRange{Int}}
end


"""
    Locus()
    Locus(position::UnitRange{Int})
    Locus(position::UnitRange{Int}, strand::Char)

"""
Locus() = Locus(1:1, '.', true, true, UnitRange{Int}[])
Locus(position::UnitRange{Int}) = Locus(position, '.', true, true, UnitRange{Int}[])
Locus(position::UnitRange{Int}, strand::Char) = Locus(position, strand, true, true, UnitRange{Int}[])
Locus(position::UnitRange{Int}, strand::Char, complete_left, complete_right) = Locus(position, strand, complete_left, complete_right, UnitRange{Int}[])

Base.convert(::Type{Locus}, x::UnitRange{Int}) = Locus(x)
function Base.convert(::Type{Locus}, x::StepRange{Int, Int})
    if x.step == -1
        return Locus(x.stop:x.start, '-')
    elseif x.step == 1
        return Locus(x.start:x.stop)
    end
    throw(DomainError(x, "`x` must have a step of 1 or -1"))
end

mutable struct Chromosome{G <: AbstractGene}
    name::String
    sequence::LongDNASeq
    header::String
    genes::Vector{G}
    genedata::DataFrame
    function Chromosome{G}() where G
        new("", dna"", "", G[], DataFrame(feature = String[], locus = Locus[]))
    end
end


Chromosome(args...) = Chromosome{Gene}(args...)


struct Gene <: AbstractGene
    index::UInt
    parent::Chromosome
end


"""
    addgene!(chr::Chromosome, feature, locus; kw...)

Add gene to `chr`. `locus` can be a `Locus`, a UnitRange, or a StepRange.
"""
function addgene!(chr::Chromosome, feature, locus; kw...)
    locus = convert(Locus, locus)
    push!(chr.genedata, vcat([feature, locus], fill(missing, size(chr.genedata, 2) - 2)))
    index = UInt32(length(chr.genes) + 1)
    gene = Gene(index, chr)
    for (k, v) in kw
        Base.setproperty!(gene, k, v)
    end
    push!(chr.genes, gene)
    return gene
end


"""
    delete!(gene::Gene)

Delete `gene` from `gene.parent`. Warning: does not work when broadcasted! Use
`delete!(::AbstractVector{Gene}) instead.`
"""
function Base.delete!(gene::Gene)
    i = gene.index
    deleterows!(gene.parent.genedata, i)
    deleteat!(gene.parent.genes, lastindex(gene.parent.genes))
    nothing
end


"""
    delete!(genes::AbstractArray{Gene, 1})

Delete all genes in `genes` from `genes[1].parent`.
"""
function Base.delete!(genes::AbstractArray{Gene, 1})
    indices = genes.index
    DataFrames.deleterows!(genes.parent.genedata, Int.(indices))
    lastindices = length(genes.parent.genes) - length(indices) + 1 : length(genes.parent.genes)
    deleteat!(genes.parent.genes, lastindices)
    nothing
end


function Base.propertynames(gene::AbstractGene)
    names(gene.parent.genedata)
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
    pushproperty!(gene::AbstractGene, name::Symbol, x::T)

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
function pushproperty!(gene::AbstractGene, name::Symbol, x::T; forceany = true) where T
    if hasproperty(gene.parent.genedata, name)
        C = eltype(gene.parent.genedata[!, name])
        if T <: C
            if ismissing(gene.parent.genedata[gene.index, name])
                gene.parent.genedata[gene.index, name] = x
            else
                gene.parent.genedata[!, name] = vectorise(gene.parent.genedata[!, name])
                push!(gene.parent.genedata[gene.index, name], x)
            end
        elseif Vector{T} <: C
            if ismissing(gene.parent.genedata[gene.index, name])
                gene.parent.genedata[gene.index, name] = [x]
            else
                push!(gene.parent.genedata[gene.index, name], x)
            end
        elseif forceany && !(C <: AbstractVector)
            if ismissing(gene.parent.genedata[gene.index, name])
                gene.parent.genedata[!, name] = convert(Vector{Any}, gene.parent.genedata[!, name])
                gene.parent.genedata[gene.index, name] = x
            else
                gene.parent.genedata[!, name] = vectorise(convert(Vector{Any}, gene.parent.genedata[! ,name]))
                push!(gene.parent.genedata[gene.index, name], x)
            end
        elseif forceany && C <: AbstractVector
            if ismissing(gene.parent.genedata[gene.index, name])
                gene.parent.genedata[!, name] = convert(Vector{Any}, gene.parent.genedata[!, name])
                gene.parent.genedata[gene.index, name] = [x]
            else
                @error "This shouldn't happen"
            end
        else
            @error "Tried to add a '$T' to '$name::$(typeof(gene.parent.genedata[!, name]))'"
        end
    else
        s = size(gene.parent.genedata, 1)
        gene.parent.genedata[!, name] = Vector{Union{Missing, T}}(missing, s)
        gene.parent.genedata[gene.index, name] = x
    end
    return x
end


"""
    get(g::Gene, key, default)

Retrieve `key` from `g`. If `key` is missing, return `default`.
"""
function Base.get(g::Gene, key, default)
    p = getproperty(g, key)
    if ismissing(p)
        return default
    else
        return p
    end
end


"""
    sequence(gene::AbstractGene)

Return genomic sequence for `gene`.
"""
function sequence(gene::AbstractGene)
    ifelse(gene.locus.strand == '-',
        reverse_complement(gene.parent.sequence[gene.parent.genedata[gene.index, :locus].position]),
        gene.parent.sequence[gene.parent.genedata[gene.index, :locus].position])
end


function Base.length(gene::AbstractGene)
    length(gene.parent.genedata[gene.index, :locus].position)
end


"""
    iscomplement(gene::Abstract)

Return `true` if `gene.locus.compliment == '-'`, otherwise return `false`.
"""
iscomplement(gene::AbstractGene) = gene.parent.genedata[gene.index, :locus].strand == '-'


function appendstring(field, v)
    s = ""
    if v isa Bool
        s *= "\n" * join(fill(' ', 21)) * "/$field"
    elseif v isa Union{Number, Symbol}
        s *= "\n" * join(fill(' ', 21)) * "/$field=" * string(v)
    else
        v = replace(v, "\n" => "\n" * join(fill(' ', 21)))
        s *= "\n" * join(fill(' ', 21)) * "/$field=\"" * v * "\""
    end
    return s
end

function Base.show(io::IO, gene::AbstractGene)
    buf = IOBuffer()
    print(buf, "     " * rpad(string(gene.feature), 16, ' '))
    print(buf, gene.locus)
    for field in names(gene.parent.genedata)
        if field in [:feature, :locus]
            continue
        end
        v = gene.parent.genedata[gene.index, field]
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


Base.show(io::IO, mime::MIME"text/plain", genes::Array{Gene, 1}) = Base.show(io, genes)
function Base.show(io::IO, genes::Array{Gene, 1})
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
        s *= "order(" * join([join((r.start, r.stop), "..") for r in locus.order], ",") * ")"
    else
        s *= join((locus.position.start, locus.position.stop), "..")
    end
    !locus.complete_right   && (s *= "<")
    locus.strand == '-'     && (s *= ")")
    print(io, s)
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


function printgbk(chrs::AbstractVector{C}) where {C <: Chromosome}
    io = IOBuffer()
    printgbk(io, chrs)
end
function printgbk(io::IO, chrs::AbstractVector{C}) where {C <: Chromosome}
    for chr in chrs
        printgbk(io, chr)
    end
    return io
end
function printgbk(chr::C) where {C <: Chromosome}
    io = IOBuffer()
    printgbk(io, chr)
end
function printgbk(io::IO, chr::C) where {C <: Chromosome}
    println(io, chr.header)
    println(io, rpad("FEATURES", 21, ' '), "Location/Qualifiers")
    println(io, chr.genes)
    println(io, "ORIGIN")
    formatsequence(chr.sequence, io)
    println(io)
    println(io, "//")
    return io
end


Base.copy(chr::Chromosome) = deepcopy(chr)
Base.copy(gene::AbstractGene) = deepcopy(gene)


function Base.isless(l1::Locus, l2::Locus)
    if (l1.position.start < l2.position.start) || (l1.position.start == l2.position.start && l1.position.stop > l2.position.stop)
        return true
    end
    return false
end
Base.isless(g1::Gene, g2::Gene) = ((g1.locus == g2.locus) && (g1.feature == "gene" && g2.feature != "gene")) || (g1.locus < g2.locus)


function Base.:(==)(x::Locus, y::Locus)
    (x.position == y.position) && (x.strand == y.strand) && (x.complete_left == y.complete_left) && (x.complete_right == y.complete_right) && (x.order == y.order)
end


function Base.sort!(G::AbstractArray{Gene}; args...)
    sort!(G[1].parent.genedata; cols = (:locus, :feature), args...)
end
