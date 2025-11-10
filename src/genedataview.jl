struct GeneDataView{G <: AbstractGene} <: AbstractArray{G, 1}
    parent::Vector{Record{G}}
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


function Base.getproperty(genes::G, name::Symbol) where {G <: AbstractVector{Gene}}
    if name in fieldnames(G)
        getfield(genes, name)
    else
        GeneDataView(parent.(genes), index.(genes), name)
    end
end


function Base.setproperty!(gene::G, name::Symbol, x::T) where {G <: AbstractGene, T}
    if hasproperty(parent(gene).genedata, name)
        parent(gene).genedata[index(gene), name] = x
    else
        s = isempty(parent(gene).genedata) ? length(parent(gene).genes) : size(parent(gene).genedata, 1)
        parent(gene).genedata[!, name] = Vector{Union{Missing, T}}(missing, s)
        parent(gene).genedata[index(gene), name] = x
    end
    return x
end

function Base.getproperty(locus::SpanLocus{T}, s::Symbol) where T
    if s == :start
        return first(locus.position)
    elseif s == :stop
        return last(locus.position)
    elseif s == :strand
        return '+'
    end
    return getfield(locus, s)
end

function Base.getproperty(locus::PointLocus{T}, s::Symbol) where T
    if s == :start
        return first(locus.position)
    elseif s == :stop
        return last(locus.position)
    elseif s == :strand
        return '+'
    end
    return getfield(locus, s)
end

function Base.getproperty(locus::Complement{T}, s::Symbol) where T
    if s == :start
        return first(locus.loc.position)
    elseif s == :stop
        return last(locus.loc.position)
    elseif s == :position
        return locus.loc.position
    elseif s == :strand
        return '-'
    end
    return getfield(locus, s)
end

function Base.getproperty(locus::Union{Join{T}, Order{T}}, s::Symbol) where T
    if s == :start
        return first(locus.loc[1].position)
    elseif s == :stop
        return last(locus.loc[end].position)
    elseif s == :position
        return Iterators.flatten(loc.position for loc in locus.loc)
    elseif s == :strand
        return locus.loc[1].strand
    end
    return getfield(locus, s)
end

Base.propertynames(::AbstractLocus) = (:start, :stop, :positition, :strand)


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
