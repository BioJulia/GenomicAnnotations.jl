# GenBank Reader

struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    io::S
end

"""
    GenBank.Reader(input::IO)

Create a data reader of the GenBank file format.
"""
function Reader(input::IO)
    if input isa TranscodingStream
        Reader(input)
    else
        stream = TranscodingStreams.NoopStream(input)
        Reader(stream)
    end
end

function Base.eltype(::Type{<:Reader})
    return Record
end

function Base.close(reader::Reader)
    close(reader.io)
end

function Base.open(::Type{<:Reader}, input::AbstractString)
    if input[end-2:end] == ".gz"
        return Reader(GzipDecompressorStream(open(input)))
    else
        return Reader(TranscodingStreams.NoopStream(open(input)))
    end
end

function BioGenerics.IO.stream(reader::Reader)
    reader.io
end

function Base.iterate(reader::Reader, nextone::Record = Record())
    if BioGenerics.IO.tryread!(reader, nextone) === nothing
        return nothing
    end
    return nextone, Record()
end

function BioGenerics.IO.tryread!(reader::Reader, output)
    record = parsechromosome!(reader.io, output)
    return record
end

"""
Return the LOCUS entry of the header.
"""
function parseheader(header::String)
    lines = split(header, "\n")
    firstline = lines[findfirst(line -> occursin("LOCUS", line), lines)]
    locus = split(firstline, r" +")[2]
    return locus
end


"""
Parse lines encoding genomic position, returning the feature as a `Symbol`, and an instance of `AbstractLocus`.
"""
function parseposition(line::String)
    f1 = findfirst(!=(' '), line)
    f2 = f1 + findfirst(==(' '), @view(line[f1:end])) - 1
    feature = Symbol(@view(line[f1:(f2-1)]))
    p1 = f2 + findfirst(!=(' '), @view(line[f2:end])) - 1
    loc = Locus(@view(line[p1:end]))
    return feature, loc
end


"""
Parse footer (sequence) portion of a GenBank file, returning a `String`.
"""
function filterseq(io::IOBuffer)
    line = String(take!(io))
    String(replace(filter(x -> !(isspace(x) || isnumeric(x)), line), r"[^ATGCNatgcn]" => "n"))
end


"""
    parsechromosome(stream::IO)

Parse and return one chromosome entry, and the line number that it ends at.
"""
function parsechromosome!(stream::IO, record::Record{G}) where G <: AbstractGene
    eof(stream) && return nothing
    iobuffer = IOBuffer()
    isheader = true
    isfooter = false
    spanning = false
    position_spanning = false
    qualifier = String("")
    content = String("")
    header = ""

    feature = :source
    loc = Locus("1..1")


    linecount = 0
    while !eof(stream)
        line = readline(stream)
        linecount += 1

        if length(line) == 0
            continue
        end

        # Catch cases where there's no header
        if linecount == 1 && @view(line[1:9]) == "     gene"
            isheader = false
        end

        ### HEADER
        if isheader && @view(line[1:8]) == "FEATURES"
            record.header = String(take!(iobuffer))
            isheader = false

        elseif isheader
            linecount == 1 ? print(iobuffer, line) : print(iobuffer, '\n', line)

        # Check if the footer has been reached
        elseif !isheader && !isfooter && (occursin(r"^BASE COUNT", line) || occursin(r"^ORIGIN", line))
            # Stop parsing the file when the list of genes is over
            isfooter = true
            iobuffer = IOBuffer()

        ### BODY
        elseif !isheader && !isfooter
            # if position_spanning && occursin(r"  /", line)
            #     position_spanning = false
            #     spanningline = filter(x -> x != ' ', String(take!(iobuffer)))
            #     try
            #         feature, locus = parseposition(spanningline)
            #     catch
            #         println(spanningline)
            #         println(line)
            #         @error "parseposition(spanningline) failed at line $linecount"
            #     end
            # elseif position_spanning
            #     print(iobuffer, line)
            # end
            if occursin(r"^ {5}\S", line)
                spanning = false
                try
                    feature, loc = parseposition(line)
                catch
                    println(line)
                    @error "parseposition(line) failed at line $linecount"
                end
                addgene!(record, feature, loc)
            elseif !spanning && occursin(r"^ +/", line)
                if occursin("=", line)
                    if occursin("=\"", line)
                        q1 = findfirst(==('/'), line) + 1
                        q2 = q1 + findfirst(==('='), @view(line[q1:end])) - 1
                        qualifier = line[q1:(q2-1)]
                        c1 = q2 + findfirst(!=('='), @view(line[q2:end]))
                        content = line[c1:end-1]
                    else
                        q1 = findfirst(==('/'), line) + 1
                        q2 = q1 + findfirst(==('='), @view(line[q1:end])) - 1
                        qualifier = line[q1:(q2-1)]
                        c1 = q2 + findfirst(!=('='), @view(line[q2:end])) - 1
                        content = line[c1:end]
                        try
                            tmpcontent = Meta.parse(content)
                            tmpcontent isa Expr && throw(Meta.ParseError)
                            content = tmpcontent
                        catch
                            content = Symbol(content)
                        end
                    end

                    if occursin(r"\".*[^\"]$", line)
                        spanning = true
                    end

                    isempty(names(record.genedata)) ?
                        (record.genedata[!, Symbol(qualifier)] = Union{Missing, typeof(content)}[content]) :
                        pushproperty!(record.genes[end], Symbol(qualifier), content)

                else
                    # Qualifiers without a value assigned to them end up here
                    q1 = findfirst(==('/'), line) + 1
                    qualifier = line[q1:end]
                    pushproperty!(record.genes[end], Symbol(qualifier), true)
                end
            elseif spanning
                if line[end] == '"'
                    content = line[22:end-1]
                    spanning = false
                else
                    content = line[22:end]
                end
                if eltype(record.genedata[!, Symbol(qualifier)]).b <: AbstractArray
                    record.genedata[!, Symbol(qualifier)][end][end] = oneline(Base.getproperty(record.genes[end], Symbol(qualifier))[end] * "\n" * content)
                else
                    Base.setproperty!(record.genes[end], Symbol(qualifier), oneline(Base.getproperty(record.genes[end], Symbol(qualifier)) * "\n" * content))
                end
            end

        ### FOOTER
        elseif isfooter
            if line == "//"
                break
            end
            if line[1] == ' '
                print(iobuffer, line)
            end
        end
    end
    if isempty(record.header) && isempty(record.genes) && isempty(record.sequence)
        return nothing
    end
    record.name = parseheader(record.header)
    record.sequence = LongDNA{4}(filterseq(iobuffer))
    return record
end


function parsegene!(record, geneindex, bytes, safe)
    newlines = findall(==(UInt8('\n')), bytes)
    fp = 4 + findfirst(==(UInt8(' ')), bytes[6:newlines[1]])
    feature = Symbol(bytes[6:fp])
    loc = Locus(String(bytes[22:(newlines[1]-1)]))
    qualifiers = Dict{Symbol, Any}()
    qualifier = :locus_tag
    for p in [(start = x[1], stop = x[2]) for x in zip(newlines[1:end-1], newlines[2:end])]
        line = bytes[p.start+1:p.stop-1]
        if line[22] == UInt8('/')
            qp = findfirst(==(UInt8('=')), line)
            if !isnothing(qp)
                qualifier = Symbol(line[23:qp-1])
                if line[qp+1] == UInt8('"')
                    ### Strings
                    if line[end] != UInt8('"')
                        content = String(line[qp+2:end])
                    else
                        content = String(line[qp+2:end-1])
                    end
                elseif line[qp+1] == UInt8('(')
                    ### Anticodon strings
                    content = String(line[qp+1:end])
                else
                    ### Ints etc.
                    content = String(line[qp:end])
                    try
                        tmpcontent = Meta.parse(content)
                        tmpcontent isa Expr && throw(Meta.ParseError)
                        content = tmpcontent
                    catch
                        content = Symbol(content)
                    end
                end
            else
                ### pseudo etc.
                qualifier = Symbol(line[23:end])
                content = true
            end
            ### add feature
            if haskey(qualifiers, qualifier)
                if typeof(qualifiers[qualifier]) isa AbstractVector
                    push!(qualifiers[qualifier], content)
                else
                    qualifiers[qualifier] = [qualifiers[qualifier], content]
                end
            else
                qualifiers[qualifier] = content
            end
        else
            ### Spanning
            if line[end] == UInt8('"')
                content = String(line[22:end-1])
            else
                content = String(line[22:end])
            end
            if haskey(qualifiers, qualifier)
                if qualifiers[qualifier] isa AbstractVector
                    qualifiers[qualifier][end] = oneline(qualifiers[qualifier][end] * "\n" * content)
                else
                    qualifiers[qualifier] = oneline(qualifiers[qualifier] * "\n" * content)
                end
            end
        end
    end
    lock(safe)
    # lock(safe) do
    try
        for (k, v) in qualifiers
            pushproperty!(record, geneindex, k, v)
        end
        # addgene!(record, feature, loc, geneindex; qualifiers...)
    finally
        unlock(safe)
    end
    Gene(record, geneindex, loc, feature)
end


function parsegene!(record, geneindex, bytes, safe)
    newlines = findall(==(UInt8('\n')), bytes)
    fp = 4 + findfirst(==(UInt8(' ')), bytes[6:newlines[1]])
    feature = Symbol(bytes[6:fp])
    loc = Locus(String(bytes[22:(newlines[1]-1)]))
    qualifier = :locus_tag
    content_str = ""
    content_int = 0
    content_symbol = :CDS
    type = :str
    for (i, p) in enumerate([(start = x[1], stop = x[2]) for x in zip(newlines[1:end-1], newlines[2:end])])
        line = bytes[p.start+1:p.stop-1]
        if line[22] == UInt8('/')
            if i != 1
                lock(safe)
                try
                    if type == :str
                        pushproperty!(record, geneindex, qualifier, content_str)
                    elseif type == :int
                        pushproperty!(record, geneindex, qualifier, content_int)
                    elseif type == :symbol
                        pushproperty!(record, geneindex, qualifier, content_symbol)
                    elseif type == :bool
                        pushproperty!(record, geneindex, qualifier, true)
                    end
                finally
                    unlock(safe)
                end
            end
            qp = findfirst(==(UInt8('=')), line)
            if !isnothing(qp)
                qualifier = Symbol(line[23:qp-1])
                if line[qp+1] == UInt8('"')
                    ### Strings
                    type = :str
                    if line[end] != UInt8('"')
                        content_str = String(line[qp+2:end])
                    else
                        content_str = String(line[qp+2:end-1])
                    end
                elseif line[qp+1] == UInt8('(')
                    type = :str
                    ### Anticodon strings
                    content_str = String(line[qp+1:end])
                else
                    ### Ints etc.
                    type = :int
                    content_str = String(line[qp:end])
                    try
                        tmpcontent = Meta.parse(content_str)
                        tmpcontent isa Expr && throw(Meta.ParseError)
                        content_int = tmpcontent
                    catch
                        type = :symbol
                        content_symbol = Symbol(content_str)
                    end
                end
            else
                ### pseudo etc.
                type = :bool
                qualifier = Symbol(line[23:end])
                # content = true
            end
        else
            ### Spanning
            if line[end] == UInt8('"')
                content_str = content_str * String(line[22:end-1])
            else
                content_str = content_str * String(line[22:end])
            end
        end
    end
    lock(safe)
    try
        if type == :str
            pushproperty!(record, geneindex, qualifier, content_str)
        elseif type == :int
            pushproperty!(record, geneindex, qualifier, content_int)
        elseif type == :symbol
            pushproperty!(record, geneindex, qualifier, content_symbol)
        elseif type == :bool
            pushproperty!(record, geneindex, qualifier, true)
        end
    finally
        unlock(safe)
    end
    Gene(record, geneindex, loc, feature)
end


function parsechromosome!(stream::IO, record::Record{G}) where G <: AbstractGene
    eof(stream) && return record
    iobuffer = IOBuffer()
    isfooter = false
    spanning = false
    position_spanning = false
    qualifier = String("")
    content = String("")
    header = ""
    missingheader = false

    feature = :source
    loc = Locus("1..1")

    record.genedata.locus_tag = Vector{Union{Missing, String}}()

    linecount = 0

    ### HEADER
    while !eof(stream)
        line = readline(stream)
        linecount += 1

        if length(line) == 0
            continue
        end

        # Catch cases where there's no header
        if linecount == 1 && @view(line[1:9]) == "     gene"
            missingheader = true
            break
        end

        if @view(line[1:8]) == "FEATURES"
            record.header = String(take!(iobuffer))
            break
        else
            linecount == 1 ? print(iobuffer, line) : print(iobuffer, '\n', line)
        end
    end
        
    ### BODY
    ### First, search the stream for the positions of all genes
    stream_gene_positions = Int[]
    while !eof(stream)
        linestartpos = position(stream)
        line = readline(stream)
        if occursin(r"^ {5}\S", line)
            ### new gene
            push!(stream_gene_positions, linestartpos)
        elseif occursin(r"^BASE COUNT", line) || occursin(r"^ORIGIN", line) || occursin(r"^CONTIG", line)
            push!(stream_gene_positions, linestartpos)
            break
        end
    end
    offset = stream_gene_positions[1]
    seek(stream, offset)


    ### BODY
    ### Next, parse each gene
    stream_gene_startstop = [(start = x[1], stop = x[2]) for x in zip(stream_gene_positions[1:end-1], stream_gene_positions[2:end])]
    safe_seek = ReentrantLock()
    safe_push = ReentrantLock()
    genebytes = Vector{UInt8}(undef, stream_gene_positions[end] - stream_gene_positions[1])
    readbytes!(stream, genebytes, length(genebytes))
    Threads.@threads for (geneindex, p) in collect(enumerate(stream_gene_startstop))
        bytes = view(genebytes, (p.start-offset+1):(p.stop-offset+1))
        parsegene!(record, geneindex, bytes, safe_push)
    end
    # sort!(record.genes)
    seek(stream, stream_gene_positions[end])

    ### FOOTER
    while !eof(stream)
        line = readline(stream)
        if line == "//"
            break
        end
        if line[1] == ' '
            print(iobuffer, line)
        end
    end
    if isempty(record.header) && isempty(record.genes) && isempty(record.sequence)
        return nothing
    end
    record.name = parseheader(record.header)
    record.sequence = LongDNA{4}(filterseq(iobuffer))
    return record
end


### Spawn
function parsechromosome!(stream::IO, record::Record{G}) where G <: AbstractGene
    eof(stream) && return nothing
    iobuffer = IOBuffer()
    isfooter = false
    spanning = false
    position_spanning = false
    qualifier = String("")
    content = String("")
    header = ""
    missingheader = false

    feature = :source
    loc = Locus("1..1")


    linecount = 0

    ### HEADER
    while !eof(stream)
        line = readline(stream)
        linecount += 1

        if length(line) == 0
            continue
        end

        # Catch cases where there's no header
        if linecount == 1 && occursin(r"^     gene", line)
            missingheader = true
            break
        end

        if occursin(r"^FEATURES", line)
            record.header = String(take!(iobuffer))
            break
        else
            linecount == 1 ? print(iobuffer, line) : print(iobuffer, '\n', line)
        end
    end
        
    ### BODY
    ### First, search the stream for the positions of all genes
    stream_gene_positions = Int[]
    # while !eof(stream)
    #     linestartpos = position(stream)
    #     line = readline(stream)
    #     if occursin(r"^ {5}\S", line)
    #         ### new gene
    #         push!(stream_gene_positions, linestartpos)
    #     elseif occursin(r"^BASE COUNT", line) || occursin(r"^ORIGIN", line) || occursin(r"^CONTIG", line)
    #         push!(stream_gene_positions, linestartpos)
    #         break
    #     end
    # end

    while !eof(stream)
        linestartpos = position(stream)
        line = readuntil(stream, UInt8('\n'), keep = true)
        if line[6] != UInt8(' ') && all(==(UInt8(' ')), line[1:5])
            push!(stream_gene_positions, linestartpos)
        elseif line[1] != UInt8(' ') && occursin(r"^(BASE COUNT|ORIGIN|CONTIG)", String(line))
            push!(stream_gene_positions, linestartpos)
            break
        end
    end

    offset = stream_gene_positions[1]
    seek(stream, offset)

    ### BODY
    ### Next, parse each gene
    stream_gene_startstop = [(start = x[1], stop = x[2]) for x in zip(stream_gene_positions[1:end-1], stream_gene_positions[2:end])]
    safe_push = ReentrantLock()
    record.genes = [Gene(record, 1, Locus("1"), :CDS) for _ in eachindex(stream_gene_startstop)]
    # bytes = [UInt8[] for _ in eachindex(record.genes)]
    # for (geneindex, p) in collect(enumerate(stream_gene_startstop))
    #     seek(stream, p.start)
    #     # readbytes!(stream, bytes[geneindex], p.stop - p.start)
    # end
    genebytes = Vector{UInt8}(undef, stream_gene_positions[end] - stream_gene_positions[1])
    readbytes!(stream, genebytes, length(genebytes))
    record.genedata.locus_tag = missings(String, length(stream_gene_startstop))
    for (geneindex, p) in collect(enumerate(stream_gene_startstop))
        bytes = genebytes[(p.start-offset+1):(p.stop-offset)]
        record.genes[geneindex] = parsegene!(record, geneindex, bytes, safe_push)
    end
    seek(stream, stream_gene_positions[end])

    ### FOOTER
    while !eof(stream)
        line = readline(stream)
        if line == "//"
            break
        end
        if line[1] == ' '
            print(iobuffer, line)
        end
    end
    if isempty(record.header) && isempty(record.genes) && isempty(record.sequence)
        return nothing
    end
    record.name = parseheader(record.header)
    record.sequence = LongDNA{4}(filterseq(iobuffer))
    return record
end


function parsegene!(record, geneindex, bytes)
    newlines = findall(==(UInt8('\n')), bytes)
    fp = 4 + findfirst(==(UInt8(' ')), bytes[6:newlines[1]])
    feature = Symbol(bytes[6:fp])
    loc = Locus(String(bytes[22:(newlines[1]-1)]))
    qualifiers = Dict{Symbol, Any}()
    qualifier = :locus_tag
    addgene!(record, feature, loc, geneindex)
    for p in [(start = x[1], stop = x[2]) for x in zip(newlines[1:end-1], newlines[2:end])]
        line = bytes[p.start+1:p.stop-1]
        if line[22] == UInt8('/')
            qp = findfirst(==(UInt8('=')), line)
            if !isnothing(qp)
                qualifier = Symbol(line[23:qp-1])
                if line[qp+1] == UInt8('"')
                    ### Strings
                    if line[end] != UInt8('"')
                        content = String(line[qp+2:end])
                    else
                        content = String(line[qp+2:end-1])
                    end
                elseif line[qp+1] == UInt8('(')
                    ### Anticodon strings
                    content = String(line[qp+1:end])
                else
                    ### Ints etc.
                    content = String(line[qp+1:end])
                    try
                        tmpcontent = Meta.parse(content)
                        tmpcontent isa Expr && throw(Meta.ParseError)
                        content = tmpcontent
                    catch
                        content = Symbol(content)
                    end
                end
            else
                ### pseudo etc.
                qualifier = Symbol(line[23:end])
                content = true
            end
            ### add feature
            if haskey(qualifiers, qualifier)
                if typeof(qualifiers[qualifier]) isa AbstractVector
                    push!(qualifiers[qualifier], content)
                else
                    qualifiers[qualifier] = [qualifiers[qualifier], content]
                end
            else
                qualifiers[qualifier] = content
            end
        else
            ### Spanning
            if line[end] == UInt8('"')
                content = String(line[22:end-1])
            else
                content = String(line[22:end])
            end
            if haskey(qualifiers, qualifier)
                if qualifiers[qualifier] isa AbstractVector
                    qualifiers[qualifier][end] = oneline(qualifiers[qualifier][end] * "\n" * content)
                else
                    qualifiers[qualifier] = oneline(qualifiers[qualifier] * "\n" * content)
                end
            end
        end
    end
    for (k, v) in qualifiers
        pushproperty!(record, geneindex, k, v)
    end
end

function parsechromosome!(stream::IO, record::Record{G}) where G <: AbstractGene
    eof(stream) && return nothing
    iobuffer = IOBuffer()
    isheader = true
    isfooter = false
    spanning = false
    position_spanning = false
    qualifier = String("")
    content = String("")
    header = ""

    feature = :source
    loc = Locus("1..1")


    linecount = 0
    ### HEADER
    while !eof(stream)
        line = readline(stream)
        linecount += 1

        if length(line) == 0
            continue
        end

        # Catch cases where there's no header
        if linecount == 1 && occursin(r"^     gene", line)
            missingheader = true
            break
        end

        if occursin(r"^FEATURES", line)
            record.header = String(take!(iobuffer))
            break
        else
            linecount == 1 ? print(iobuffer, line) : print(iobuffer, '\n', line)
        end
    end

    ### BODY
    genebuffer = IOBuffer()
    geneindex = 0
    while !eof(stream)
        line = readline(stream)
        linecount += 1

        if length(line) == 0
            continue
        end

        if occursin(r"^     \S", line)
            if geneindex > 0
                parsegene!(record, geneindex, take!(genebuffer))
            end
            geneindex += 1
            println(genebuffer, line)
        elseif occursin(r"^(BASE COUNT|ORIGIN)", line)
            break
        elseif line[1] != ' '
            nothing
        else
            println(genebuffer, line)
        end
    end

    ### FOOTER
    iobuffer = IOBuffer()
    while !eof(stream)
        line = readline(stream)
        if line == "//"
            break
        end
        print(iobuffer, line)
    end
    close(stream)
    record.name = parseheader(record.header)
    record.sequence = LongDNA{4}(filterseq(iobuffer))
    return record
end
