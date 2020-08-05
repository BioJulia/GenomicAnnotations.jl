# I/O

## Input
Annotation files are read with `readgbk(input)` or `readgff(input)`. Currently these assume that the file follows either standard [GenBank format](http://www.insdc.org/files/feature_table.html#7.1.2), or [GFF3](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/). If a file can't be read it is likely not conforming to the specifications (but feel free to submit a bug report).

```@docs
readgbk
readgff
```

## Output
Annotations can be printed with GenBank formatting using `printgbk`, and as GFF3 with `printgff`. In the REPL, instances of `Gene` are displayed as they would be in the annotation file.

```@docs
printgbk
printgff
```
