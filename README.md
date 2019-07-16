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


## Usage
GenBank files are read with `readgbk(gbkfile)`, which returns a vector of `Chromosome`s. `gbkfile` can be gzipped as long as the filename ends in ".gz". If we're only interested in the first chromosome in `example.gbk` we only need to store the first element.
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

Accessing properties that haven't been stored will return missing. For this reason, it often makes more sense to use `get()` than to access the property directly.
```julia
# chr.genes[2].pseudo returns missing, so this will throw an error
if chr.genes[2].pseudo
    println("Gene 2 is a pseudogene")
end

# ... but this works:
if get(chr.genes[2], :pseudo, false)
    println("Gene 2 is a pseudogene")
end
```

The macro `@genes` can be used to filter through the annotations. The keyword `gene` is used to refer to the individual `Gene`s. `@genes` can also be used to modify annotations.
```julia
@genes(chr, :feature == "CDS")  # Returns all coding regions
@genes(chr, iscds)              # Shorthand for ':feature == "CDS"'
@genes(chr, length(gene) > 300) # Returns all features longer than 300 nt
@genes(chr, iscomplement(gene)) # Returns all features on the complement strand

# The following two expressions are equivalent:
@genes(chr, :feature == "CDS", length(gene) > 300)
@genes(chr, (:feature == "CDS") && (length(gene) > 300))

@genes(chr, :locus_tag == "tag03")[1].pseudo = true
delete!(@genes(chr, :pseudo)) # Delete all psudogenes. This works for all genes
                              # where `:pseudo` is `true`, and ignores genes
                              # where it is `false` or `missing`
delete!(@genes(chr, length(gene) <= 60)) # Delete all genes 60 nt or shorter

# Some short-hand forms are available to make life easier:
#     `iscds` expands to `:feature == "CDS"`, and
#     `get(s::Symbol, default)` expands to `get(gene, s, default)`
# The following two are thus equivalent:
@genes(chr, :feature == "CDS", occursin("glycoprotein", get(gene, :product, "")))
@genes(chr, iscds, occursin("glycoprotein", get(:product, "")))
```

Gene sequences can be accessed with `sequence(gene)`. For example, the following code will write the translated sequences of all protein-coding genes to a file:
```julia
using BioSequences
writer = FASTA.Writer(open("proteins.fasta", "w"))
for gene in @genes(chr, :feature == "CDS")
    aaseq = translate(sequence(gene))
    write(writer, FASTA.record(gene.locus_tag, get(gene, :product, ""), aaseq))
end
close(writer)
```

Genes can be added using `addgene!`, and `sort!` can be used to make sure that the resulting annotations are in the correct order.
```julia
newgene = addgene!(chr, "regulatory", 670:677)
newgene.locus_tag = "reg02"
sort!(chr.genes)
```

After modifying the annotations, `printgbk(io, chr)` can be used to write them to a file.
```julia
open("updated.gbk", "w") do f
    printgbk(f, chr)
end
```
