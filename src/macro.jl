function genes_helper!(ex)
    if ex isa Expr
        if ex.head == :call
            if ex.args[1] != :get
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
accessed in the expression as `gene`. Accessing properties of the returned list
of genes returns a view, which can be altered.

Some short-hand forms are available to make life easier:
    `iscds` expands to `:feature == "CDS"`, and
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
@genes(chr, :feature == "CDS", length(gene) > 300)
@genes(chr, (:feature == "CDS") && (length(gene) > 300))
```

`@genes` returns a `Vector{Gene}`. Attributes can be accessed with dot-syntax, and can be assigned to
```julia
@genes(chr, :locus_tag == "tag03")[1].pseudo = true
@genes(chr, iscds, ismissing(:gene)).gene .= "unknown"
```
"""
macro genes(chr, exs...)
    exs = Any[ex for ex in exs]
    for (i, ex) in enumerate(exs)
        # Special cases
        if ex == :iscds
            exs[i] = :(gene.feature == "CDS")
        elseif ex == :(!iscds)
            exs[i] = :(gene.feature != "CDS")
        # Expressions
        elseif ex isa Expr
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
        $hits = falses(size($chr.genes, 1))
        for ($i, gene) in enumerate($chr.genes)
            local h = $ex
            $hits[$i] = ismissing(h) ? false : h
        end
        $chr.genes[$hits]
    end)
end
