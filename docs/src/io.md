# I/O

## Input
Annotation files are read with `GenBank.Reader` and `GFF.Reader`. Currently these assume that the file follows either standard [GenBank format](http://www.insdc.org/files/feature_table.html#7.1.2), or [GFF3](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/). Any metadata in GFF3 files, apart from the header, is ignored.
```julia
open(GenBank.Reader, "example.gbk") do reader
    for record in reader
        do_something()
    end
end
```
`readgbk(input)` and `readgff(input)` are aliases for `collect(open(GenBank.Reader, input))` and `collect(open(GFF.Reader, input))`, respectively.

```@docs
GenBank.Reader
GFF.Reader
```

## Output
Annotations can be printed with GenBank formatting using `GenBank.Writer`, and as GFF3 with `GFF.Writer`. Headers are not automatically converted between formats; `GFF.Writer` only prints the header of the first `Record`, and only if it starts with a `#`, while `GenBank.Writer` prints a default header if the stored one starts with `#`.

```@docs
GenBank.Writer
GFF.Writer
```

In the REPL, instances of `Gene` are displayed as they would be in the annotation file.
