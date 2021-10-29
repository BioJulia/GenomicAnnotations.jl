"""
    oneline(s)

Remove newlines, replacing them with spaces when it seems appropriate.
"""
function oneline(v::AbstractString)
    R = split(v, '\n')
    buf = IOBuffer()
    for (i, r) in enumerate(R)
        if i != lastindex(R) && (r[end] == '-' || occursin(' ', r))
            print(buf, r, " ")
        else
            print(buf, r)
        end
    end
    String(take!(buf))
end
oneline(v::Any) = v


function findbreak(v, t)
    i = t == 59 ?
        findlast('\n', v) :
        1
    j = i+t-1
    j >= lastindex(v) && return nothing
    x = findlast(c -> c == ' ' || c == '-', v[1:j])
    if isnothing(x)
        return j:j+1
    elseif x == ' '
        return x-1:x+1
    else
        return x:x+1
    end
end


"""
    multiline(v, s::Symbol)

Return a `String` with newlines and spaces added so that it conforms to the GenBank line width.
"""
function multiline(v, s)
    v = string(v)
    s = string(s)
    if length(v) + length(s) > 55 && !occursin('\n', v)
        b = findbreak(v, 55 - length(s))
        v = v[1:b.start] * "\n" * v[b.stop:end]
        b = findbreak(v, 59)
        while !isnothing(b) && b.start < length(v)
            v = v[1:b.start] * "\n" * v[b.stop:end]
            b = findbreak(v, 59)
        end
    end
    return v
end


"""
    readgbk(input)

Parse GenBank-formatted file, returning a `Vector{GenBank.Record}`. File names ending in ".gz" are assumed to be gzipped and are decompressed.
"""
function readgbk(input)
    collect(open(GenBank.Reader, input))
end


"""
    readgff(input)

Parse GFF3-formatted file, returning a `Vector{Record}`. File names ending in ".gz" are assumed to be gzipped and are decompressed.
"""
function readgff(input)
    collect(open(GFF.Reader, input))
end


function relocate_gene(gene, pos, chrlen; reverse = false)
    oldpos = locus(gene).position
    newposition = if reverse
            (chrlen - mod1(oldpos.stop - pos - 1, chrlen)) : (chrlen - mod1(oldpos.start - pos - 1, chrlen))
        elseif length(oldpos) < chrlen
            mod1(oldpos.start - pos + 1, chrlen) : mod1(oldpos.stop - pos + 1, chrlen)
        else
            oldpos
        end
    newstrand = if !reverse
            locus(gene).strand
        elseif locus(gene).strand == '+'
            '-'
        elseif locus(gene).strand == '-'
            '+'
        else
            locus(gene).strand
        end
    order = locus(gene).order
    join = locus(gene).join
    newlocus = Locus(newposition, newstrand, locus(gene).complete_right, locus(gene).complete_left, order, join)
    Gene(parent(gene), index(gene), newlocus, feature(gene))
end


"""
    reorder!(chr, pos = 1; reverse = false)

Reorder `chr` so that `pos` becomes the first position. If `reverse` is `true`, the reverse genome is reversed.
"""
function reorder(chr, pos = 1; reverse = false)
    newchr = deepcopy(chr)
    reorder!(newchr, pos; reverse = reverse)
end
function reorder!(chr, pos = 1; reverse = false)
    ## Reverse sequence
    seq = reverse ?
        BioSequences.reverse_complement(chr.sequence[pos+1:end] * chr.sequence[1:pos]) :
        chr.sequence[pos:end] * chr.sequence[1:pos-1]
    ## Recalculate gene loci
    chrlen = length(chr.sequence)
    genes = Gene[]
    for gene in chr.genes
        newgene = relocate_gene(gene, pos, chrlen; reverse = reverse)
        push!(genes, newgene)
    end
    chr.genes = genes
    sort!(chr.genes)
    chr.sequence = seq
    chr
end