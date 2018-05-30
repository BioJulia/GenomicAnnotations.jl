# __precompile__(true)
module GenomicAnnotations

using DataFrames

export Chromosome, Gene, AbstractGene
export readgbk, genesequence, iscomplement, addgene!, printgbk
export @genes

include("types.jl")
include("readgbk.jl")
include("macro.jl")

end #module
