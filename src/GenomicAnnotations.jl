module GenomicAnnotations

using TranscodingStreams
using BioGenerics
using DataFrames
using BioSequences

export GenBank, GFF, EMBL
export Gene, AbstractGene, GeneDataView
export sequence, iscomplement, iscomplete, ismultilocus, addgene!, pushproperty!
export feature, index, locus, locus!, position

export AbstractLocus
export SpanLocus, ClosedSpan, OpenSpan, OpenRightSpan, OpenLeftSpan
export PointLocus, SingleNucleotide, BetweenNucleotides
export Join, Order, Complement
export Locus

export relative_position

export @genes, upstream, downstream, neighbours
export readgbk, readgff
export reorder, reorder!

include("record.jl")
include("genedataview.jl")
include("macro.jl")
include("utils.jl")


include("GenBank/GenBank.jl")
import .GenBank
include("GFF/GFF.jl")
import .GFF
include("EMBL/EMBL.jl")
import .EMBL

end #module
