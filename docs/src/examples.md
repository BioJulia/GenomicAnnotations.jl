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

open(GenBank.Writer, "updated_genome.gbk") do w
    for chr in chrs
        write(w, chr)
    end
end
```


## Converting between formats
Annotations can be read from one file format and written as another. Converting between the supported human-readable formats (GenBank and EMBL) or the tab-delimited formats (GFF3 and GTF) will likely work out of the box, but converting from a human-readable format to a tab-delimited format, or vice versa, may need some human intervention. Currently, GenomicAnnotations does not make any attempt to rename columns or perform any sanity checks to ensure that the resulting file meets specifications. Refer to the respective format specifications for details on what attributes need to be included, etc. Notably, GenBank and GFF3 headers do not contain the same information, and GTF files lack a header altogether, thus all information in the header is lost when saving annotations as another format. GTF files also do not allow the inclusion of sequence data, unlike GFF3.
Below is a simple example script that demonstrates some of the changes that need to be made. If your use-case includes more complex features, such as multi-exon genes, you will likely need to make more changes as part of the convertion.
```julia
using GenomicAnnotations
using DataFrames
chrs = readgbk("genome.gbk")
open(GFF.Writer, "genome.gff") do w
    for chr in chrs
        # GenBank features often contain features that are not usually included
        # in GFF3 files, so let's remove some:
        cols_to_remove = intersect(["translation", "mol_type", "organism"], names(chr.genedata))
        if !isempty(cols_to_remove)
            chr.genedata = chr.genedata[:, Not(cols_to_remove)]
        end
        # The GenBank format uses the :source feature to store metadata about
        # the record, but in GFF3, :source is the name of the column which
        # contains the sequence name. Thus, we need to change the GenBank
        # :source to the GFF3 equivalent :region. According to the GenBank
        # specification, :source is mandatory, but it's best to be safe and
        # check that it's really there:
        source_entries = @genes(chr, source)
        if !isempty(source_entries)
            for source in source_entries
                region = feature!(source, :region)
                region.Name = chr.name
                region.ID = chr.name
                # GenBank files include information about circularity of a
                # contig in its header, but in the GFF3 format this
                # information is encoded in the "Is_circular" attribute of
                # the first :region feature:
                if occursin("circular", chr.header)
                    region.Is_circular = true
                end
            end
        else
            # If the :source feature is missing, we'll have to create a :region
            # from scratch:
            addgene!(chr, "region", 1:length(chr.sequence);
                Name = chr.name,
                ID = chr.name,
                Is_circular = occursin("circular", chr.header))
            sort!(chr.genes)
        end
        # Most features, such as :CDS or :tRNA, have a corresponding :gene
        # that it belongs to. In GFF3, this hierarchical relationship is shown
        # using the "ID" and "Parent" attributes. Here, we set the "ID"
        # attribute of all :gene features to match their "locus_tag", and then
        # set the "Parent" attributes of all non-:gene features to match the
        # "ID" of their respective :gene, if there is one:
        gene_features = @genes(chr, gene)
        gene_features.ID .= gene_features.locus_tag
        for gene in @genes(chr, !gene)
            if get(gene, :locus_tag, "missing") in skipmissing(gene_features.locus_tag)
                gene.Parent = gene.locus_tag
            end
        end
        write(w, chr)
    end
end
```
