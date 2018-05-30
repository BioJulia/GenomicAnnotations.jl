function genes_helper!(ex)
    if ex isa Expr
        if ex.head == :call
            for i in eachindex(ex.args)[2:end]
                if ex.args[i] isa QuoteNode
                    ex.args[i] = :(Base.getproperty(gene, $(ex.args[i])))
                else
                    genes_helper!(ex.args[i])
                end
            end
        else
            for i in eachindex(ex.args)
                if ex.args[i] isa QuoteNode
                    ex.args[i] = :(Base.getproperty(gene, $(ex.args[i])))
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
accessed in the expression as `gene` (see example below). Missing values are
treated as `false`.
```jldoctest
julia> chromosome = readgbk("example.gbk")
Chromosome 'example' (5028 bp) with 6 annotations

julia> @genes(chromosome, :feature == "CDS") |> length
3

julia> @genes(chromosome, length(gene) < 500)
     CDS             3..206
                     /db_xref="GI:1293614"
                     /locus_tag="tag01"
                     /codon_start="3"
                     /product="TCP1-beta"
                     /protein_id="AAA98665.1"
```

Accessing properties of the returned list of genes returns a view, which can be
altered.
```jldoctest
julia> @genes(chromosome, ismissing(:gene)) |> length
2

julia> @genes(chromosome, ismissing(:gene)).gene .= "Unknown";

julia> @genes(chromosome, ismissing(:gene)) |> length
0
```
"""
macro genes(chr, exs...)
    exs = Any[ex for ex in exs]
    for (i, ex) in enumerate(exs)
        if ex isa Expr
            genes_helper!(ex)
        elseif ex isa QuoteNode
            exs[i] = :(Base.getproperty(gene, $ex))
        end
    end
    ex = Expr(:&&, exs...)
    hits = gensym()
    i = gensym()
    return esc(quote
        $hits = falses(size($chr.genes, 1))
        for ($i, gene) in enumerate($chr.genes)
            local h = $ex
            $hits[$i] = ismissing(h) ? false : h
        end
        $chr.genes[$hits]
    end)
end
