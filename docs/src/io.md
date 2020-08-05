# I/O

## Input
Annotation files are read with `readgbk(input)` or `readgff(input)`. Currently these assume that the file follows either standard [GenBank format](http://www.insdc.org/files/feature_table.html#7.1.2), or [GFF3](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/). Any metadata in GFF3 files, apart from the header, is ignored.

```@docs
readgbk
readgff
```

## Output
Annotations can be printed with GenBank formatting using `printgbk`, and as GFF3 with `printgff`. Headers are not automatically converted between formats; `printgff` only prints the header of the first `Chromosome`, and only if it starts with a `#`, while `printgbk` prints a default header if the stored one starts with `#`.

```@docs
printgbk
printgff
```

In the REPL, instances of `Gene` are displayed as they would be in the annotation file.
