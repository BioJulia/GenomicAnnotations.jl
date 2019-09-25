# I/O

## Input
Annotation files are read with `readgbk(pathtofile)`. Currently this assumes that the file follows standard [GenBank format](http://www.insdc.org/files/feature_table.html#7.1.2).

```@docs
readgbk(file)
```

## Output
Annotations can be printed as with GenBank formatting using `printgbk`. In the REPL, instances of `Gene` are displayed as they would be in the annotation file.

```@docs
printgbk
```
