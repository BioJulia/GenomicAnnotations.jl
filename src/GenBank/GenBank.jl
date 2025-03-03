module GenBank

using TranscodingStreams
using BioSequences
using BioGenerics
using DataFrames
using CodecZlib

using ..GenomicAnnotations
import ..GenomicAnnotations: Record, Gene, AbstractGene, GeneDataView, Locus, GenBankFormat
import ..GenomicAnnotations: sequence, iscomplement, iscomplete, addgene!, pushproperty!, feature, index, locus, oneline
export sequence, iscomplement, iscomplete, feature, index, locus

include("reader.jl")
include("writer.jl")

end #module
