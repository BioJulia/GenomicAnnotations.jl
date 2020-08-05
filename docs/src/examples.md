# Examples
## Adding chromosome name to all locus tags
When iterating over genes, the parent chromosome can be accessed with `parent(::Gene)`.
```julia
using GenomicAnnotations
chrs = readgbk("genome.gbk")
for gene in @genes(chrs)
    gene.locus_tag = string(parent(gene).name, "_", gene.locus_tag)
end
printgbk("updated_genome.gbk", chrs)
```
## Adding qualifiers
GenomicAnnotations supports arbitrary qualifiers, so you can add any kind of information. The following script reads and adds the output from [Phobius](http://phobius.sbc.su.se/) (a predictor for transmembrane helices) to the annotations.
```julia
using GenomicAnnotations
chrs = readgbk("genome.gbk")

function addphobius!(chr, file)
    @progress for line in readlines(file)
        m = match(r"^(\w+) +(\d+) +", line)
        if m != nothing
            locus_tag = m[1]
            tmds = parse(Int, m[2])
            @genes(chr, CDS, :locus_tag == locus_tag).phobius .= tmds
        end
    end
end

addphobius!(chrs, "phobius.txt")

printgbk("updated_genome.gbk", chrs)
```


## Converting between formats
Note that GenBank and GFF3 headers do not contain the same information, thus all information in the header is lost when saving annotations as another format.
```julia
using GenomicAnnotations
chrs = readgbk("genome.gbk")
printgff("genome.gff", chrs)
```
