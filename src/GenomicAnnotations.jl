__precompile__(true)
module GenomicAnnotations

using DataFrames
using BioSequences
using GZip

export Chromosome, Gene, AbstractGene, GeneDataView, Locus
export readgbk, sequence, iscomplement, iscomplete, addgene!, pushproperty!, printgbk
export feature, index, locus
export @genes, upstream, downstream, neighbours
export readgff, printgff

include("types.jl")
include("genedataview.jl")
include("readgbk.jl")
include("macro.jl")

end #module
