struct GeneDataView{G <: AbstractGene} <: AbstractArray{G, 1}
    parent::Chromosome{G}
    indices::Vector{UInt}
    property::Symbol
end


function Base.getproperty(gene::G, name::Symbol) where {G <: AbstractGene}
    if name in fieldnames(G)
        return getfield(gene, name)
    elseif name in propertynames(gene)
        return gene.parent.genedata[gene.index, name]
    else
        return missing
    end
end


function Base.getproperty(gv::GeneDataView{G}, name::Symbol) where {G <: AbstractGene}
    return getfield(gv, name)
end


function Base.getproperty(genes::AbstractArray{G, 1}, name::Symbol) where {G <: AbstractGene}
    # @assert all(getfield.(genes, :parent) .== Ref(getfield(genes[1], :parent))) "Not all members of `genes` come from the same parent"
    if name == :index
        return [getfield(g, :index) for g in genes]
    elseif name == :parent
        return getfield(genes[1], name)
    else
        return GeneDataView(genes[1].parent, genes.index, name)
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


### AbstractArray functions
Base.size(gv::GeneDataView) = size(gv.indices)
Base.getindex(gv::GeneDataView, i::Int) = getproperty(gv.parent.genes[gv.indices[i]], gv.property)
function Base.getindex(gv::GeneDataView, I::AbstractArray)
    if haskey(gv.parent.genedata, gv.property)
        return gv.parent.genedata[gv.indices[I], gv.property]
    else
        println(1)
        return fill(missing, length(gv))
    end
end
Base.setindex!(gv::GeneDataView, v, i::Int) = setindex!(gv.parent.genedata[gv.property], v, i)
Base.setindex!(gv::GeneDataView, v, I::Vararg{Int, N}) where N = setindex!(gv.parent.genedata[gv.property], v, i)
