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

The `Locus` of a `Gene` retrieved with `locus(gene)`. The `Locus` itself is immutable, but can be updated with `locus!(gene, newlocus)`. For simplicity, `position(gene)` is shorthand for `locus(gene).position`.
```julia
# Create a new Locus, copying all fields of the old one but shifting the position by 1
oldloc = locus(gene)
locus!(gene, Locus(oldloc.position .+ 1, oldloc.strand, oldloc.complete_left, oldloc.complete_right, oldloc.order, oldloc.join))

# Access the genomic positions of all genes
position.(chr.genes)
```

The macro `@genes` can be used to filter through the annotations (see [`@genes`](@ref)). The keyword `gene` is used to refer to the individual `Gene`s. `@genes` can also be used to modify annotations.
```julia
@genes(chr, length(gene) > 300) # Returns all features longer than 300 nt
```

Gene sequences can be accessed with `sequence(gene)`. For example, the following code will write the translated sequences of all protein-coding genes in `chr` to a file:
```julia
using BioSequences
using FASTX
open(FASTA.Writer, "proteins.fasta") do w
    for gene in @genes(chr, CDS)
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
