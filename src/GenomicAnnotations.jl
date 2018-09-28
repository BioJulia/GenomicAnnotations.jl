# __precompile__(true)
module GenomicAnnotations

using DataFrames

export Chromosome, Gene, AbstractGene, GeneDataView, Locus
export readgbk, genesequence, iscomplement, addgene!, printgbk
export @genes

include("types.jl")
include("genedataview.jl")
include("readgbk.jl")
include("macro.jl")

end #module
