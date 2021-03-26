# GenomicAnnotations.jl
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://BioJulia.github.io/GenomicAnnotations.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://BioJulia.github.io/GenomicAnnotations.jl/stable)

## Description
GenomicAnnotations is a package for reading, modifying, and writing genomic annotations in the GenBank and GFF3 file formats.

## Installation
```julia
julia>]
pkg> add GenomicAnnotations
```


## Usage
GenBank and GFF3 files are read with `readgbk(input)` and `readgff(input)`, which return vectors of `Record`s. `input` can be an `IOStream` or a file path. GZipped data  is unzipped automatically if a filename ending in ".gz" is passed as `input`. If we're only interested in the first chromosome in `example.gbk` we only need to store the first element.
```julia
chr = readgbk("test/example.gbk")[1]
```

`Record`s have five fields, `name`, `header`, `genes`, `genedata`, and `sequence`. The `name` is read from the `header`, which is stored as a string. The annotation data is stored in `genedata`, but generally you should use `genes` to access that data. For example, it can be used to iterate over annotations, and to modify them.
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

The macro `@genes` can be used to filter through the annotations. The macro takes a `Record` or a `Vector{Record}`, followed by any number of expressions that will be evaluated for each gene. The keyword `gene` is used to refer to the individual `Gene`s. `@genes` can also be used to modify annotations. Gene attributes can be referred to using `Symbol`s.
```julia
@genes(chr, feature(gene) == "CDS")  # Returns all coding regions
@genes(chr, length(gene) > 300) # Returns all features longer than 300 nt
@genes(chr, iscomplement(gene)) # Returns all features on the complement strand
@genes(chr, ismissing(:product)) # Returns all features for which the attribute "product" has not been set

# Some short-hand forms are available to make life easier:
#     `CDS` expands to `feature(gene) == "CDS"`, and
#     `get(s::Symbol, default)` expands to `get(gene, s, default)`
# The following two are thus equivalent:
@genes(chr, feature(gene) == "CDS", occursin("glycoprotein", get(gene, :product, "")))
@genes(chr,                   CDS,  occursin("glycoprotein", get(      :product, "")))

# All arguments have to evaluate to `true` for a gene to be included, so the following expressions are equivalent:
@genes(chr, feature(gene) == "CDS", length(gene) > 300)
@genes(chr, (feature(gene) == "CDS") && (length(gene) > 300))

# `@genes` returns a `Vector{Gene}`. Attributes can be accessed with dot-syntax, and can be assigned to:
@genes(chr, :locus_tag == "tag03")[1].pseudo = true
@genes(chr, CDS, ismissing(:gene)).gene .= "unknown"
```

Gene sequences can be accessed with `sequence(gene)`. For example, the following code will write the translated sequences of all protein-coding genes to a file:
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

Individual genes, and `Vector{Gene}`s are printed in GBK format. To include the GBK header and the nucleotide sequence, `printgbk(io, chr)` can be used to write them to a file. `printgff(io, chr)` prints the annotations as GFF3, in which case the GenBank header is lost.
```julia
println(chr.genes[1])
println(@genes(chr, CDS))

open(GenBank.Writer, "updated.gbk") do w
    write(w, chr)
end
```
