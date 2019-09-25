# GenomicAnnotations.jl

## Description
GenomicAnnotations is a package for reading, modifying, and writing genomic annotations in the GenBank format.


## Installation
```julia
julia>]
pkg> add GenomicAnnotations
```
or
```julia
using Pkg
Pkg.add("GenomicAnnotations")
```


## Examples
GenBank files are read with `readgbk(pathtofile)`, which returns a vector of `Chromosome`s. `gbkfile` can be gzipped as long as the filename ends in ".gz". If we're only interested in the first chromosome in `example.gbk` we only need to store the first element.
```julia
chr = readgbk("test/example.gbk")[1]
```

`Chromosome`s have five fields, `name`, `header`, `genes`, `genedata`, and `sequence`. The `name` is read from the `header`, which is stored as a string. The annotation data is stored in `genedata`, but generally you should use `genes` to access that data. For example, it can be used to iterate over annotations, and to modify them.
```julia
for gene in chr.genes
    gene.locus_tag = "$(chr.name)_$(gene.locus_tag)"
end

chr.genes[2].locus_tag = "test123"
```

The macro `@genes` can be used to filter through the annotations. The keyword `gene` is used to refer to the individual `Gene`s. `@genes` can also be used to modify annotations.
```julia
@genes(chr, length(gene) > 300) # Returns all features longer than 300 nt
```

Gene sequences can be accessed with `sequence(gene)`. For example, the following code will write the translated sequences of all protein-coding genes to a file:
```julia
using FASTX
writer = FASTA.Writer(open("proteins.fasta", "w"))
for gene in @genes(chr, iscds)
    aaseq = translate(sequence(gene))
    write(writer, FASTA.record(gene.locus_tag, get(gene, :product, ""), aaseq))
end
close(writer)
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

Individual genes, and `Vector{Gene}`s are printed in GBK format. To include the GBK header and the nucleotide sequence, `printgbk(io, chr)` can be used to write them to a file.
```julia
println(chr.genes[1])
println(@genes(chr, iscds))

open("updated.gbk", "w") do f
    printgbk(f, chr)
end
```
