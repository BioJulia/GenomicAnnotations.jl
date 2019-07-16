struct GeneDataView{G <: AbstractGene} <: AbstractArray{G, 1}
    parent::Chromosome{G}
    indices::Vector{UInt}
    property::Symbol
end


function Base.getproperty(gene::G, name::Symbol) where {G <: AbstractGene}
    if name in (:parent, :index)
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
    if hasproperty(gene.parent.genedata, name)
        gene.parent.genedata[gene.index, name] = x
    else
        s = size(gene.parent.genedata, 1)
        gene.parent.genedata[!, name] = Vector{Union{Missing, T}}(missing, s)
        gene.parent.genedata[gene.index, name] = x
    end
    return x
end


Base.size(gv::GeneDataView) = size(gv.indices)
Base.getindex(gv::GeneDataView, i::Int) = getproperty(gv.parent.genes[gv.indices[i]], gv.property)
Base.getindex(gv::GeneDataView, I::AbstractArray) = GeneDataView(gv.parent, gv.indices[I], gv.property)
Base.setindex!(gv::GeneDataView, v, i::Int) = (Base.setproperty!(gv.parent.genes[i], gv.property, v))
Base.similar(gv::GeneDataView{Gene}) = GeneDataView(gv.parent, gv.indices, gv.property)
Base.copy(gv::GeneDataView{Gene}) = GeneDataView(gv.parent, gv.indices, gv.property)
Base.view(gv::GeneDataView{Gene}, I) = GeneDataView(gv.parent, gv.indices[I], gv.property)


function Base.fill!(gv::GeneDataView{Gene}, x)
    if hasproperty(gv.parent.genedata, gv.property)
        xT = convert(eltype(gv.parent.genedata[!, gv.property]), x)
    else
        xT = x
    end
    for I in gv.indices
        @inbounds gv[I] = x
    end
    gv
end
