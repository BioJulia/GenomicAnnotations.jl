# __precompile__(true)
module GenomicAnnotations

using DataFrames
using BioSequences

export Chromosome, Gene, AbstractGene, GeneDataView, Locus
export readgbk, genesequence, iscomplement, addgene!, pushproperty!, printgbk
export @genes

include("types.jl")
include("genedataview.jl")
include("readgbk.jl")
include("macro.jl")

end #module
