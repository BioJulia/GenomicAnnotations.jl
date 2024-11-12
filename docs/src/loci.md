# [Representing genomci loci](@id Loci)

The easiest way to create a locus is to use the constructor `Locus(s)`, which takes an `AbstractString` `s` and parses it as a GenBank locus string as defined here: https://www.insdc.org/submitting-standards/feature-table/#3.4. Note that remote entry descriptors have not been implemented.

## Internal representation
Since v0.4.0, genomic loci are represented using instances of `AbstractLocus`. Simple descriptors are represented with `PointLocus{T}` and `SpanLocus{T}`, where `T` is an `AbstractDescriptor`:

| GenBank string | GenomicAnnotations representation | Description |
| --- | --- | --- |
| 1   | `PointLocus{SingleNucleotide}(1)` | Refers to a single nucleotide. |
| 1^2 | `PointLocus{BetweenNucleotides}(1)` | Refers to the internucleotide space immediately after position 1. |
| 10..20 | `SpanLocus{ClosedSpan}(10:20)` | Denotes a closed sequence span. |
| 10..>20 | `SpanLocus{OpenRightSpan}(10:20) | Denotes a sequence span where the right side is open, i.e. the end-point is undefined but earliest at position 20. |
| <10..20 | `SpanLocus{OpenLeftSpan}(10:20) | The left end-point is undefined. |
| <10..>20 | `SpanLocus{OpenSpan}(10:20) | Both end-points are undefined. |

These can be wrapped in `Complement` for loci on the complement strand, e.g. `Complement(SpanLocus{ClosedSpan}(10:20))` representing "complement(10..20)". Simplified constructors are provided for all `AbstractDescriptor`s, e.g. `ClosedSpan(1:10) == SpanLocus(1:10, ClosedSpan)`.

Compound loci are represented with `Join` and `Order`. Both types have a single field, `loc` which contains any number of simple descriptors. They can be wrapped with complement, as can the individual elements in `loc`.

```julia
Locus("complement(join(10..20,30..>40))") isa Complement{Join{SpanLocus{ClosedSpan}, SpanLocus{OpenRightSpan}}}
```