# Filtering: the @genes macro

A useful tool provided by GenomicAnnotations is the macro `@genes`. It is used to filter through annotations, for example to look at only at coding sequences or rRNAs, which can then be modified or iterated over:
```julia
# Print locus tags of all coding sequences longer than 1000 nt, that are not pseudo genes
for gene in @genes(chr, iscds, length(gene) > 1000, ! :pseudo)
    println(gene.locus_tag)
end
```

```@docs
@genes
```
