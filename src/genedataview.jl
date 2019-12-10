struct GeneDataView{G <: AbstractGene} <: AbstractArray{G, 1}
    parent::Vector{Chromosome{G}}
    indices::Vector{UInt}
    property::Symbol
end


function Base.getproperty(gene::G, name::Symbol) where {G <: AbstractGene}
    if name in propertynames(gene)
        return parent(gene).genedata[index(gene), name]
    else
        return missing
    end
end


function Base.getproperty(gv::GeneDataView{G}, name::Symbol) where {G <: AbstractGene}
    return getfield(gv, name)
end


function Base.getproperty(genes::AbstractArray{G, 1}, name::Symbol) where {G <: AbstractGene}
    GeneDataView(parent.(genes), getfield.(genes, Ref(:index)), name)
end


function Base.setproperty!(gene::G, name::Symbol, x::T) where {G <: AbstractGene, T}
    if hasproperty(parent(gene).genedata, name)
        parent(gene).genedata[index(gene), name] = x
    else
        s = size(parent(gene).genedata, 1)
        parent(gene).genedata[!, name] = Vector{Union{Missing, T}}(missing, s)
        parent(gene).genedata[index(gene), name] = x
    end
    return x
end


Base.size(gv::GeneDataView) = size(gv.indices)
Base.getindex(gv::GeneDataView, i::Int) = getproperty(gv.parent[i].genes[gv.indices[i]], gv.property)
Base.getindex(gv::GeneDataView, I::AbstractArray) = GeneDataView(gv.parent[I], gv.indices[I], gv.property)
Base.setindex!(gv::GeneDataView, v, i::Int) = (Base.setproperty!(gv.parent[i].genes[gv.indices[i]], gv.property, v))
Base.similar(gv::GeneDataView{G}) where {G <: AbstractGene} = GeneDataView(gv.parent, gv.indices, gv.property)
Base.copy(gv::GeneDataView{G}) where {G <: AbstractGene} = GeneDataView(gv.parent, gv.indices, gv.property)
Base.view(gv::GeneDataView{G}, I) where {G <: AbstractGene} = GeneDataView(gv.parent[I], gv.indices[I], gv.property)


function Base.fill!(gv::GeneDataView{G}, x) where {G <: AbstractGene}
    chrs = unique(gv.parent)
    for chr in chrs
        if hasproperty(chr.genedata, gv.property)
            xT = convert(eltype(chr.genedata[!, gv.property]), x)
        else
            xT = x
        end
    end
    for i in eachindex(gv.indices)
        @inbounds gv[i] = x
    end
    gv
end
