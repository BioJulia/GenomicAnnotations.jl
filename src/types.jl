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
    excluding::Vector{UnitRange{Int}}
end


Locus() = Locus(1:1, '.', true, true, UnitRange{Int}[])


mutable struct Chromosome{G <: AbstractGene}
    name::String
    sequence
    header::String
    genes::Vector{G}
    genedata::DataFrame
    function Chromosome{G}() where G
        new("", "", "", G[], DataFrame(feature = String[], locus = Locus[]))
    end
end


Chromosome(args...) = Chromosome{Gene}(args...)


struct Gene <: AbstractGene
    index::UInt
    parent::Chromosome
end


"""
    addgene!(chr::Chromosome, locus::Locus; kw...)

Add gene to `chr`.
"""
function addgene!(chr::Chromosome, feature, locus::Locus; kw...)
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

Delete all genes in `genes` from `gene[1].parent`.
"""
function Base.delete!(genes::AbstractArray{Gene, 1})
    indices = genes.index
    deleterows!(genes[1].parent.genedata, indices)
    lastindices = length(gene.parent.genes) - length(indices) + 1 : length(gene.parent.genes)
    deleteat!(gene.parent.genes, lastindices)
    nothing
end


#=
function Base.getproperty(gene::G, name::Symbol) where {G <: AbstractGene}
    if name in fieldnames(G)
        return getfield(gene, name)
    elseif haskey(gene.parent.genedata, name)
        return gene.parent.genedata[gene.index, name]
    end
    return missing
end


function Base.getproperty(genes::AbstractArray{G, 1}, name::Symbol) where {G <: AbstractGene}
    if name in fieldnames(G)
        return getfield.(genes, name)
    elseif haskey(genes[1].parent.genedata, name)
        return view(genes[1].parent.genedata[name], [gene.index for gene in genes])
    else
        return fill(missing, length(genes))
    end
end


function Base.setproperty!(gene::G, name::Symbol, x::T) where {G <: AbstractGene, T}
    if haskey(gene.parent.genedata, name)
        gene.parent.genedata[gene.index, name] = x
    else
        s = size(gene.parent.genedata, 1)
        gene.parent.genedata[name] = Vector{Union{Missing, T}}(missing, s)
        gene.parent.genedata[gene.index, name] = x
    end
    return x
end
=#


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

Add a property to `gene`, similarly to `Base.setproperty!(::gene)`, but
transform the property to store a vector instead of overwriting existing data.
"""
function pushproperty!(gene::AbstractGene, name::Symbol, x::T) where T
    if haskey(gene.parent.genedata, name)
        C = eltype(gene.parent.genedata[name])
        if T <: C
            if ismissing(gene.parent.genedata[gene.index, name])
                gene.parent.genedata[gene.index, name] = x
            else
                gene.parent.genedata[name] = vectorise(gene.parent.genedata[name])
                push!(gene.parent.genedata[gene.index, name], x)
            end
        elseif Vector{T} <: C
            if ismissing(gene.parent.genedata[gene.index, name])
                gene.parent.genedata[gene.index, name] = [x]
            else
                push!(gene.parent.genedata[gene.index, name], x)
            end
        else
            @error "Tried to add a '$T' to a '$(typeof(gene.parent.genedata[name]))'"
        end
    else
        s = size(gene.parent.genedata, 1)
        gene.parent.genedata[name] = Vector{Union{Missing, T}}(missing, s)
        gene.parent.genedata[gene.index, name] = x
    end
    return x
end


"""
    genesequence(gene::AbstractGene)

Return genomic sequence for `gene`.
"""
function genesequence(gene::AbstractGene)
    gene.parent.sequence[gene.parent.genedata[gene.index, :position]]
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
    s = "     " * rpad(string(gene.feature), 16, ' ')
    iscomplement(gene) && (s *= "complement(")
    !gene.locus.complete_left && (s *= "<")
    s *= "$(gene.locus.position.start)..$(gene.locus.position.stop)"
    !gene.locus.complete_right && (s *= ">")
    iscomplement(gene) && (s *= ")")
    for field in names(gene.parent.genedata)
        if field in [:feature, :locus]
            continue
        end
        v = gene.parent.genedata[gene.index, field]
        if !ismissing(v)
            if v isa AbstractVector
                for i in eachindex(v)
                    s *= appendstring(field, v[i])
                end
            else
                s *= appendstring(field, v)
            end
        end
    end
    println(io, s)
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
    if !locus.complete_left
        s = s * ">"
    end
    s = s * string(locus.position.start)
    if length(locus.excluding) > 0
        for gap in locus.excluding
            s = s * ".."
        end
    end
    if length(locus.position) > 1
        s = s * ".." * string(locus.position.stop)
    end
    if !locus.complete_right
        s = s * "<"
    end
    if locus.strand == '-'
        s = "complement(" * s * ")"
    end
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
    if (l1.position.start < l2.position.start) || (l1.position.start == l2.position.start && l1.position.stop >= l2.position.stop)
        return true
    end
    return false
end
Base.isless(g1::Gene, g2::Gene) = g1.locus < g2.locus


function Base.sort!(G::AbstractArray{Gene}; args...)
    sort!(G[1].parent.genedata; cols = (:locus, :feature), args...)
end
