"""
    oneline(s)

Remove newlines, replacing them with spaces when it seems appropriate.
"""
function oneline(v::AbstractString)
    R = split(v, '\n')
    buf = IOBuffer()
    for (i, r) in enumerate(R)
        if i != lastindex(R) && (r[end] == '-' || occursin(' ', r))
            print(buf, strip(r), " ")
        else
            print(buf, r)
        end
    end
    String(take!(buf))
end
oneline(v::Any) = v


function findbreak(v, t; width = 80, margin = 21)
    i = t == width - margin + 1 ?
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
    multiline(v, s::Symbol; width = 80, margin = 21, first_line_margin = 4)

Return a `String` with newlines and spaces added so that it conforms to the specified line width, with the margin subtracted. For example, the line width for the GenBank format is 80, with a margin of 21, whereas for the EMBL format the margin is 19.
"""
function multiline(v, s; width = 80, margin = 21, first_line_margin = 3)
    v = string(v)
    s = string(s)
    if length(v) + length(s) > width - margin - first_line_margin + 1 && !occursin('\n', v)
        b = findbreak(v, width - margin - first_line_margin + 1 - length(s))
        v = v[1:b.start] * "\n" * v[b.stop:end]
        b = findbreak(v, width - margin + 1; width, margin)
        while !isnothing(b) && b.start < length(v)
            v = v[1:b.start] * "\n" * v[b.stop:end]
            b = findbreak(v, width - margin + 1; width, margin)
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
    newposition = if length(oldpos) == chrlen
            oldpos
        elseif reverse
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

function locus!(gene::AbstractGene, loc)
    chr = parent(gene)
    newgene = Gene(chr, index(gene), loc, feature(gene))
    chr.genes[index(gene)] = newgene
end


function Base.sort!(genes::Vector{Gene}; kwargs...)
    I = sortperm(genes; kwargs...)
    oldgenes = deepcopy(genes)
    for (i, gene) in enumerate(oldgenes[I])
        genes[i] = Gene(parent(gene), UInt(i), locus(gene), feature(gene))
    end
    parent(genes[1]).genedata = parent(genes[1]).genedata[I, :]
end