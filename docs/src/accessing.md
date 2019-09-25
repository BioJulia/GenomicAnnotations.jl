# Accessing and modifying annotations

# Feature
Features (genes) can be added using `addgene!`. A feature must have a feature name and a locus (position), and can have any number of additional qualifiers associated with it (see next section).
```@docs
addgene!
```

After adding a new feature, `sort!` can be used to make sure that the annotations are stored (and printed) in the order in which they occur on the chromosome:
```julia
sort!(chr)
```

Existing features can be removed using `delete!`:
```@docs
delete!(::Gene)
delete!(::AbstractVector{Gene})
```

# Qualifiers
Features can have multiple qualifiers, which can be modified using Julia's property syntax:
```julia
# Remove newspace from gene product descriptions
for gene in @genes(chr, iscds)
    replace!(gene.product, '\n' => ' ')
end
```

Properties also work on views of genes, typically generated using `@genes`:
```julia
interestinggenes = readlines("/path/to/list/of/interesting/genes.txt")
@genes(chr, iscds, :locus_tag in interestinggenes).interesting .= true
```

Sometimes features have multiple instances of the same qualifier, such genes having several EC-numbers. Assigning qualifiers with property syntax overwrites any data that was previously stored for that feature, and trying to assign a vector of values to a qualifier that is currently storing scalars will result in an error, so to safely assign qualifiers that might have more instances one can use `pushproperty!`:
```@docs
pushproperty!
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

# Sequences
The sequence of a `Chromosome` `chr` is stored in `chr.sequence`. Sequences of individual features can be read with `sequence`:
```@docs
sequence(::Gene)
```
