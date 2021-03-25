module GenomicAnnotations

using TranscodingStreams
using BioGenerics
using DataFrames
using BioSequences

export GenBank, GFF
export Gene, AbstractGene, GeneDataView, Locus
export sequence, iscomplement, iscomplete, addgene!, pushproperty!
export feature, index, locus
export @genes, upstream, downstream, neighbours
export readgbk, readgff

include("record.jl")
include("genedataview.jl")
include("macro.jl")
include("utils.jl")


include("GenBank/GenBank.jl")
import .GenBank
include("GFF/GFF.jl")
import .GFF

end #module
