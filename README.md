# GenomicAnnotations.jl

## Description
GenomicAnnotations is a package for reading, modifying, and writing genomic annotations in the GenBank format.

## Usage
GenBank files are read with `readgbk(gbkfile)`. `readgbk(gbkfile)` returns a vector of `Chromosome`s, but since `example.gbk` only contains one we only need to store the first element.
```
chr = readgbk("test/example.gbk")[1]
```

`Chromosome`s have four fields, `name`, `genes`, `genedata`, and `sequence`. The annotation data is stored in `genedata`, but generally you should use `genes` to access that data. For example, it can be used to iterate over annotations, and to modify them.
```
for gene in chr.genes
    gene.locus_tag = "$(chr.name)_$(gene.locus_tag)"
end

chr.genes[2].locus_tag = "test123"
```

Accessing properties that haven't been stored will return missing. For this reason, it often makes more sense to use `get()` than to access the property directly.
```
if chr.genes[2].pseudo
    println("Gene 2 is a pseudogene")
end
# chr.genes[2].pseudo returns missing, so this will throw an error

if get(chr.genes[2], :pseudo, false)
    println("Gene 2 is a pseudogene")
end
```

The macro `@genes` can be used to filter through the annotations. The keyword `gene` is used to refer to the individual `Gene`s. `@genes` can also be used to modify annotations.
```
@genes(chr, :feature == "CDS")
@genes(chr, length(gene) > 300)

# The following two expressions are equivalent:
@genes(chr, :feature == "CDS", length(gene) > 300)
@genes(chr, (:feature == "CDS") && (length(gene) > 300))

@genes(chr, :locus_tag == "tag03")[1].pseudo = true
delete!(@genes(chr, :pseudo))
```

After modifying the annotations, `printgbk(io, chr)` can be used to write them to a file.
```
open("updated.gbk", "w") do f
    printgbk(f, chr)
end
```
