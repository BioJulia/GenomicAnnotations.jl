# GenomicAnnotations.jl

## Description
GenomicAnnotations is a package for reading, modifying, and writing genomic annotations in the GenBank and GFF3 file formats.


## Installation
```julia
julia>]
pkg> add GenomicAnnotations
```


## Examples
GenBank and GFF3 files are read with `readgbk(input)` and `readgff(input)`, which return vectors of `Record`s. `input` can be an `IOStream` or a file path. GZipped data can be read by setting the keyword `gunzip` to true, which is done automatically if a filename ending in ".gz" is passed as `input`. If we're only interested in the first chromosome in `example.gbk` we only need to store the first record.
```julia
chr = readgbk("test/example.gbk")[1]
```
Another way to read files is to use the corresponding `Reader` directly:
```julia
open(GenBank.Reader, "test/example.gbk") do reader
    for record in reader
        println(record.name)
    end
end
```

`Record`s have five fields, `name`, `header`, `genes`, `genedata`, and `sequence`. The `name` is read from the `header`, which is stored as a string. The annotation data is stored in `genedata`, but generally you should use `genes` to access that data. For example, it can be used to iterate over annotations, and to modify them.
```julia
for gene in chr.genes
    gene.locus_tag = "$(chr.name)_$(gene.locus_tag)"
end

chr.genes[2].locus_tag = "test123"
```

The locus of a `Gene` is represented by an `AbstractLocus` (see [Loci](@ref)), which can be retrieved with `locus(gene)`. The locus of a gene can be updated with `locus!(gene, newlocus)`. The easiest way to create a locus is to use the constructor `Locus(s)`, which takes an `AbstractString` `s` and parses it as a GenBank locus string as defined here: https://www.insdc.org/submitting-standards/feature-table/#3.4. Note that remote entry descriptors have not been implemented.
```julia
# Creating a new locus
newlocus = Locus("complement(join(1..100,200..>300))")

# Assigning a new locus to a gene
locus!(gene, newlocus)
# which is equivalent to
locus!(gene, "complement(join(1..100,200..>300))")
```

For simplicity, `position(gene)` is shorthand for `locus(gene).position`. `locus(gene).position` gives an iteratable object that generates each individual position in the defined order. Thus:
```julia
loc = Locus("join(4..6,1..3)")
collect(loc.position) # Returns [4,5,6,1,2,3]
```

The macro `@genes` can be used to filter through the annotations (see [`@genes`](@ref)). The keyword `gene` is used to refer to the individual `Gene`s. `@genes` can also be used to modify annotations.
```julia
@genes(chr, length(gene) > 300) # Returns all features longer than 300 nt

@genes(chr, CDS, ismissing(:product)) .= "hypothetical product"
```

Gene sequences can be accessed with `sequence(gene)`, which returns the nucleotide sequence. If the `translate` keyword is set to `true`, the translated amino acid sequence is returned instead. By default the first codon is translated to methionine also for alternate start codons, but this behaviour can be toggled by setting `preserve_alternate_start` to `false`. No checks are made to ensure that the gene points to a valid open reading frame, so this should be done by the user. The following example will write the translated sequences of all protein-coding genes in `chr` to a file:
```julia
using BioSequences
using FASTX
open(FASTA.Writer, "proteins.fasta") do w
    for gene in @genes(chr, CDS, iscomplete(gene))
        aaseq = GenomicAnnotations.sequence(gene; translate = true)
        write(w, FASTA.Record(gene.locus_tag, get(:product, ""), aaseq))
    end
end
```

Genes can be added using `addgene!`, and `sort!` can be used to make sure that the resulting annotations are in the correct order for printing. `delete!` is used to remove genes.
```julia
newgene = addgene!(chr, "regulatory", 670:677)
newgene.locus_tag = "reg02"
sort!(chr.genes)

# Genes can be deleted. This works for all genes where `:pseudo` is `true`, and ignores genes where it is `false` or `missing`
delete!(@genes(chr, :pseudo))
# Delete all genes 60 nt or shorter
delete!(@genes(chr, length(gene) <= 60))
```

Individual genes, and `Vector{Gene}`s are printed in GBK format. To include the GBK header and the nucleotide sequence, `write(::GenBank.Writer, chr)` can be used to write them to a file. Use `GFF.Writer` instead to print the annotations as GFF3, in which case the GenBank header is lost.
```julia
println(chr.genes[1])
println(@genes(chr, CDS))

open(GenBank.Writer, "updated.gbk") do w
    write(w, chr)
end
```
