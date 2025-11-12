module GenomicAnnotations

using TranscodingStreams
using BioGenerics
using DataFrames
using BioSequences
using FASTX

export GenBank, GFF, GTF, EMBL
export Gene, AbstractGene, GeneDataView
export sequence, iscomplement, iscomplete, iscompound, addgene!, pushproperty!
export feature, feature!, index, locus, locus!, position, attributes, genedata

export AbstractLocus
export SpanLocus, ClosedSpan, OpenSpan, OpenRightSpan, OpenLeftSpan
export PointLocus, SingleNucleotide, BetweenNucleotides
export Join, Order, Complement
export Locus
export eachposition
export shift, shift!

export relative_position

export @genes, upstream, downstream, neighbours
export readgbk, readgff, readgtf, readembl
export reorder, reorder!

include("locus.jl")
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
include("GTF/GTF.jl")
import .GTF

end #module
