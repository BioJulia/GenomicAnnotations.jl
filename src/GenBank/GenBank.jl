module GenBank

using TranscodingStreams
using BioSequences
using BioGenerics
using DataFrames
using CodecZlib

import ..GenomicAnnotations: Record, Gene, AbstractGene, GeneDataView, Locus
import ..GenomicAnnotations: sequence, iscomplement, iscomplete, addgene!, pushproperty!, feature, index, locus, oneline
export sequence, iscomplement, iscomplete, feature, index, locus

include("reader.jl")
include("writer.jl")
include("header.jl")
export iscircular!

end #module
