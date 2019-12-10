function genes_helper!(ex)
    if ex isa Expr
        if ex.args[1] != :Ref
            skipgene = ex.args[1] ∉ [:(||), :(&&)]
            skipfirst = ex.head == :call
            genes_special_cases!(ex.args; skipgene = skipgene, skipfirst = skipfirst)
        end
        if ex.head == :call
            if ex.args[1] == :Ref
                # Remove Ref(), but ignore the content
                ex2 = copy(ex)
                ex.head = :ref
                empty!(ex.args)
                push!(ex.args, ex2)
            elseif ex.args[1] != :get
                for i in eachindex(ex.args)[2:end]
                    if ex.args[i] isa QuoteNode
                        ex.args[i] = :(Base.getproperty(gene, $(ex.args[i])))
                    else
                        genes_helper!(ex.args[i])
                    end
                end
            else
                # Expand `get(s::Symbol, default)` to `get(gene, s::Symbol, default)`
                if length(ex.args) == 3
                    ex.args = vcat(ex.args[1], :gene, ex.args[2:3])
                end
            end
        else
            for i in eachindex(ex.args)
                if ex.args[i] isa QuoteNode
                    if ex.head != :(.)
                        ex.args[i] = :(Base.getproperty(gene, $(ex.args[i])))
                    end
                else
                    genes_helper!(ex.args[i])
                end
            end
        end
    end
end


"""
    @genes(chr, exs...)

Iterate over and evaluate expressions in `exs` for all genes in `chr.genes`,
returning genes where all expressions evaluate to `true`. Any given symbol `s`
in the expression will be substituted for `gene.s`. The gene itself can be
accessed in the expression as `gene`. Accessing properties of the returned list
of genes returns a view, which can be altered.

Symbols and expressions escaped with `Ref()` will be ignored.

Some short-hand forms are available to make life easier:
    `CDS`, `rRNA`, and `tRNA` expand to `feature(gene) == "..."`,
    `get(s::Symbol, default)` expands to `get(gene, s, default)`

# Examples
```julia
julia> chromosome = readgbk("example.gbk")
Chromosome 'example' (5028 bp) with 6 annotations

julia> @genes(chromosome, iscds) |> length
3

julia> @genes(chromosome, length(gene) < 500)
     CDS             3..206
                     /db_xref="GI:1293614"
                     /locus_tag="tag01"
                     /codon_start="3"
                     /product="TCP1-beta"
                     /protein_id="AAA98665.1"

julia> @genes(chromosome, ismissing(:gene)) |> length
2

julia> @genes(chromosome, ismissing(:gene)).gene .= "Unknown";

julia> @genes(chromosome, ismissing(:gene)) |> length
0
```

All arguments have to evaluate to `true` for a gene to be included, so the following expressions are equivalent:
```julia
@genes(chr, feature(gene) == Ref(:CDS), length(gene) > 300)
@genes(chr, (feature(gene) == Ref(:CDS)) && (length(gene) > 300))
```

`@genes` returns a `Vector{Gene}`. Attributes can be accessed with dot-syntax, and can be assigned to
```julia
@genes(chr, :locus_tag == "tag03")[1].pseudo = true
@genes(chr, CDS, ismissing(:gene)).gene .= "unknown"
```
"""
macro genes(chr, exs...)
    exs = Any[ex for ex in exs]
    genes_special_cases!(exs; skipgene = false, skipfirst = false)
    for (i, ex) in enumerate(exs)
        # Expressions
        if ex isa Expr
            exs[i] = :((x -> ismissing(x) ? false : x)($(exs[i])))
            genes_helper!(exs[i])
        elseif ex isa QuoteNode
            exs[i] = :((x -> ismissing(x) ? false : x)(Base.getproperty(gene, $ex)))
        end
    end
    ex = Expr(:&&, exs...)
    hits = gensym()
    i = gensym()
    return esc(quote
        if $chr isa AbstractVector
            vcat([begin
                    $hits = falses(size(c.genes, 1))
                    for ($i, gene) in enumerate(c.genes)
                        local h = $ex
                        $hits[$i] = ismissing(h) ? false : h
                    end
                    c.genes[$hits]
                end for c in $chr]...)
        else
            $hits = falses(size($chr.genes, 1))
            for ($i, gene) in enumerate($chr.genes)
                local h = $ex
                $hits[$i] = ismissing(h) ? false : h
            end
            $chr.genes[$hits]
        end
    end)
end



function genes_special_cases!(exs; skipgene = false, skipfirst = true)
    for (i, ex) in enumerate(exs)
        skipfirst && i == 1 && continue
        if ex == :CDS
            exs[i] = :(GenomicAnnotations.feature(gene) == Ref(:CDS))
        elseif ex == :(!CDS)
            exs[i] = :(GenomicAnnotations.feature(gene) != Ref(:CDS))
        elseif ex == :rRNA
            exs[i] = :(GenomicAnnotations.feature(gene) == Ref(:rRNA))
        elseif ex == :(!rRNA)
            exs[i] = :(GenomicAnnotations.feature(gene) != Ref(:rRNA))
        elseif ex == :tRNA
            exs[i] = :(GenomicAnnotations.feature(gene) == Ref(:tRNA))
        elseif ex == :(!tRNA)
            exs[i] = :(GenomicAnnotations.feature(gene) != Ref(:tRNA))
        elseif ex == :gene && !skipgene
            exs[i] = :(GenomicAnnotations.feature(gene) == Ref(:gene))
        elseif ex == :(!gene)
            exs[i] = :(GenomicAnnotations.feature(gene) != Ref(:gene))
        elseif ex == :regulatory
            exs[i] = :(GenomicAnnotations.feature(gene) == Ref(:regulatory))
        elseif ex == :(!regulatory)
            exs[i] = :(GenomicAnnotations.feature(gene) != Ref(:regulatory))
        elseif ex == :source
            exs[i] = :(GenomicAnnotations.feature(gene) == Ref(:source))
        elseif ex == :(!source)
            exs[i] = :(GenomicAnnotations.feature(gene) != Ref(:source))
        end
    end
end


# Features from the GenBank specifications without shorthand forms:
#   assembly_gap
#   C_region
#   centromere
#   D-loop
#   D_segment
#   gap
#   iDNA
#   intron
#   J_segment
#   mat_peptide
#   misc_binding
#   misc_difference
#   misc_feature
#   misc_recomb
#   misc_RNA
#   misc_structure
#   mobile_element
#   modified_base
#   mRNA
#   ncRNA
#   N_region
#   old_sequence
#   operon
#   oriT
#   polyA_site
#   precursor_RNA
#   prim_transcript
#   primer_bind
#   propeptide
#   protein_bind
#   repeat_region
#   rep_origin
#   S_region
#   sig_peptide
#   stem_loop
#   STS
#   telomere
#   tmRNA
#   transit_peptide
#   unsure
#   V_region
#   V_segment
#   variation
#   3'UTR
#   5'UTR
