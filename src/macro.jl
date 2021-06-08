function genes_helper!(ex, gene)
    if ex isa Expr && !isempty(ex.args) && ex.head != :($)
        if ex.args[1] ∉ [:upstream, :downstream, :neighbours]
            skipgene = ex.args[1] ∉ [:(||), :(&&)]
            skipfirst = ex.head == :call
            genes_special_cases!(ex.args, gene; skipgene = skipgene, skipfirst = skipfirst)
        end
        if ex.head == :call
            if ex.args[1] == :get
                if length(ex.args) == 3
                    # Expand `get(s::Symbol, default)` to `get(gene, s::Symbol, default)`
                    genes_helper!(ex.args[3], gene)
                    ex.args = vcat(ex.args[1], gene, ex.args[2:3])
                else
                    if ex.args[2] == :gene
                        ex.args[2] = gene
                    end
                    genes_helper!(ex.args[2], gene)
                    genes_helper!(ex.args[4], gene)
                end
            elseif ex.args[1] in [:upstream, :downstream, :neighbours]
                # Expand `upstream(i, f)` to `upstream(gene, f, i)`, etc.
                if length(ex.args) == 3
                    feature_function!(ex.args)
                    ex.args = vcat(ex.args[1], gene, ex.args[2:3])
                end
            else
                for i in eachindex(ex.args)[2:end]
                    if ex.args[i] isa QuoteNode
                        ex.args[i] = :(Base.getproperty($gene, $(ex.args[i])))
                    elseif ex.args[i] == :gene
                        ex.args[i] = gene
                    else
                        genes_helper!(ex.args[i], gene)
                    end
                end
            end
        elseif ex.head == :(.) && ex.args[1] == :get
            if length(ex.args[2].args) == 2
                genes_helper!(ex.args[2].args[2], gene)
                ex.args[2].args = vcat(:($gene), ex.args[2].args)
            else
                genes_helper!(ex.args[2].args[1], gene)
                genes_helper!(ex.args[2].args[3], gene)
            end
        else
            for i in eachindex(ex.args)
                if ex.args[i] isa QuoteNode
                    if ex.head != :(.)
                        ex.args[i] = :(Base.getproperty($gene, $(ex.args[i])))
                    end
                else
                    genes_helper!(ex.args[i], gene)
                end
            end
        end
    elseif ex isa Expr && ex.head === :($)
        # Allow circumventing special syntax with "$"
        ex2 = copy(ex)
        ex.head = :call
        empty!(ex.args)
        push!(ex.args, :identity)
        append!(ex.args, ex2.args)
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
    `CDS`, `rRNA`, and `tRNA` expand to `feature(gene) == "..."`,
    `get(s::Symbol, default)` expands to `get(gene, s, default)`

# Examples
```julia
julia> chromosome = readgbk("example.gbk")
Chromosome 'example' (5028 bp) with 6 annotations

julia> @genes(chromosome, CDS) |> length
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
@genes(chr, CDS, length(gene) > 300)
@genes(chr, CDS && (length(gene) > 300))
```

`@genes` returns a `Vector{Gene}`. Attributes can be accessed with dot-syntax, and can be assigned to
```julia
@genes(chr, :locus_tag == "tag03")[1].pseudo = true
@genes(chr, CDS, ismissing(:gene)).gene .= "unknown"
```

Symbols and expressions escaped with `\$` will be ignored.
```julia
d = Dict(:category1 => ["tag01", "tag02"], :category2 => ["tag03"])
@genes(chr, :locus_tag in d[\$:category1])

gene = chr.genes[5]
@genes(chr, gene == \$gene)
```
"""
macro genes(chr, exs...)
    hits = gensym()
    i = gensym()
    gene = gensym()
    exs = Any[ex for ex in exs]
    genes_special_cases!(exs, gene; skipgene = false, skipfirst = false)
    for (i, ex) in enumerate(exs)
        # Expressions
        if ex isa Expr
            exs[i] = :((x -> ismissing(x) ? false : x)($(exs[i])))
            genes_helper!(exs[i], gene)
        elseif ex isa QuoteNode
            exs[i] = :((x -> ismissing(x) ? false : x)(Base.getproperty($gene, $ex)))
        end
    end
    ex = Expr(:&&, exs...)
    return esc(quote
        if $chr isa AbstractVector{Gene}
            $hits = falses(size($chr, 1))
            for ($i, $gene) in enumerate($chr)
                local h = $ex
                $hits[$i] = ismissing(h) ? false : h
            end
            $chr[$hits]
        elseif $chr isa AbstractVector
            vcat([begin
                    $hits = falses(size(c.genes, 1))
                    for ($i, $gene) in enumerate(c.genes)
                        local h = $ex
                        $hits[$i] = ismissing(h) ? false : h
                    end
                    c.genes[$hits]
                end for c in $chr]...)
        else
            $hits = falses(size($chr.genes, 1))
            for ($i, $gene) in enumerate($chr.genes)
                local h = $ex
                $hits[$i] = ismissing(h) ? false : h
            end
            $chr.genes[$hits]
        end
    end)
end


upstream(f::Function, gene, cond::Function, i) = f(upstream(gene, cond, i))
upstream(gene, cond::Function, i::Int) = upstream(gene, cond, i:i)[1]
function upstream(gene, cond::Function, i)
    ngenes = length(parent(gene).genes)
    j = 0
    genes = Gene[]
    grange = iscomplement(gene) ?
        parent(gene).genes[mod1.(index(gene) + 1 : index(gene) + ngenes - 2, ngenes)] :
        parent(gene).genes[reverse(mod1.((index(gene) + 1) : (index(gene) + ngenes - 2), ngenes))]
    for g in grange
        if cond(g)
            j += 1
            j in i && push!(genes, g)
        end
    end
    return genes
end

downstream(f::Function, gene, cond::Function, i) = f(downstream(gene, cond, i))
downstream(gene, cond::Function, i::Int) = downstream(gene, cond, i:i)[1]
function downstream(gene, cond::Function, i)
    ngenes = length(parent(gene).genes)
    j = 0
    genes = Gene[]
    grange = !iscomplement(gene) ?
        parent(gene).genes[mod1.(index(gene) + 1 : index(gene) + ngenes - 2, ngenes)] :
        parent(gene).genes[reverse(mod1.((index(gene) + 1) : (index(gene) + ngenes - 2), ngenes))]
    for g in grange
        if cond(g)
            j += 1
            j in i && push!(genes, g)
        end
    end
    return genes
end

neighbours(f::Function, gene, cond::Function, i) = f(neighbours(gene, cond, i))
function neighbours(gene, cond::Function, i)
    vcat(upstream(gene, cond, i), downstream(gene, cond, i))
end


function feature_function!(exs; skipgene = false, skipfirst = true)
    i = 2
    ex = exs[i]
    if ex isa Symbol
        if ex == :CDS
            exs[i] = :(GenomicAnnotations.feature(g) == :CDS)
        elseif ex == :(!CDS)
            exs[i] = :(GenomicAnnotations.feature(g) != :CDS)
        elseif ex == :rRNA
            exs[i] = :(GenomicAnnotations.feature(g) == :rRNA)
        elseif ex == :(!rRNA)
            exs[i] = :(GenomicAnnotations.feature(g) != :rRNA)
        elseif ex == :tRNA
            exs[i] = :(GenomicAnnotations.feature(g) == :tRNA)
        elseif ex == :(!tRNA)
            exs[i] = :(GenomicAnnotations.feature(g) != :tRNA)
        elseif ex == :gene && !skipgene
            exs[i] = :(GenomicAnnotations.feature(g) == :gene)
        elseif ex == :(!gene)
            exs[i] = :(GenomicAnnotations.feature(g) != :gene)
        elseif ex == :regulatory
            exs[i] = :(GenomicAnnotations.feature(g) == :regulatory)
        elseif ex == :(!regulatory)
            exs[i] = :(GenomicAnnotations.feature(g) != :regulatory)
        elseif ex == :source
            exs[i] = :(GenomicAnnotations.feature(g) == :source)
        elseif ex == :(!source)
            exs[i] = :(GenomicAnnotations.feature(g) != :source)
        end
        exs[i] = Expr(:(->), :g, Expr(:block, exs[i]))
    end
end


function genes_special_cases!(exs, gene; skipgene = false, skipfirst = true)
    features = [:CDS, :rRNA, :tRNA, :mRNA, :exon, :intron, :regulatory, :source, :repeat_region]
    for (i, ex) in enumerate(exs)
        skipfirst && i == 1 && continue
        if ex isa Symbol
            if ex == :gene && !skipgene
                exs[i] = :(GenomicAnnotations.feature($gene) == $(Expr(:($), QuoteNode(:gene))))
            elseif ex in features
                exs[i] = :(GenomicAnnotations.feature($gene) == $(Expr(:($), QuoteNode(ex))))
            end
        elseif ex isa Expr && ex.head == :call && ex.args[1] == :(!)
            feature = ex.args[2]
            if feature == :gene
                exs[i] = :(GenomicAnnotations.feature($gene) != $(Expr(:($), QuoteNode(:gene))))
            elseif feature in features
                exs[i] = :(GenomicAnnotations.feature($gene) != $(Expr(:($), QuoteNode(feature))))
            end
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
