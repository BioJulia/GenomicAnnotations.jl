using GenomicAnnotations
using BioSequences
using Test

@testset "GenomicAnnotations" begin
    @testset "readgbk" begin
        s = "     gene            1..1"
        @test GenomicAnnotations.parseposition(s) == ("gene", Locus(1:1, '+'))
        s = "     gene            complement(order(3300..4037,4047..4052))"
        @test GenomicAnnotations.parseposition(s) == ("gene", Locus(3300:4052, '-', true, true, UnitRange{Int}[3300:4037, 4047:4052]))


        chrs = readgbk("example.gbk")
        @test length(chrs) == 2
        @test chrs[2].name == "plasmid1"
    end

    chr = readgbk("example.gbk")[1]

    @testset "Extended methods" begin
        @test length(chr.genes[1]) == length(sequence(chr.genes[1]))
    end

    @testset "Gene properties" begin
        @test length(propertynames(chr.genes[1])) == 14
        @test chr.genes[2].locus_tag == "tag01"
        @test (chr.genes[2].locus_tag = "tag01") == "tag01"
        @test begin
            chr.genes[2].test = "New column"
            chr.genes[2].test == "New column" && all(ismissing, chr.genes[[1,3,4,5,6]].test)
        end
        @test get(chr.genes[1], :locus_tag, "") == ""
        @test get(chr.genes[2], :locus_tag, "") == "tag01"
        @test get(chr.genes[2], :locustag, "") == ""
        @test begin
            GenomicAnnotations.pushproperty!(chr.genes[2], :db_xref, "GI:123")
            chr.genes[2].db_xref == ["GI:1293614", "GI:123"]
        end
        @test GenomicAnnotations.vectorise(Union{Missing,Int}[1,1,1]) == [[1],[1],[1]]
    end

    @testset "Iteration" begin
        @test length([g.locus_tag for g in chr.genes]) == 7
    end

    @testset "Adding/removing genes" begin
        addgene!(chr, "CDS", Locus(300:390, '+'), locus_tag = "tag04")
        @test chr.genes[end].locus_tag == "tag04"
        delete!(chr.genes[end])
        @test chr.genes[end].locus_tag == "reg01"
    end

    @testset "@genes" begin
        @test @genes(chr, :feature == "CDS") == chr.genes[[2,4,6]]
        @test @genes(chr, iscomplement(gene)) == chr.genes[[5,6,7]]
        @test @genes(chr, :feature == "CDS", !iscomplement(gene)) == chr.genes[[2,4]]
        @test @genes(chr, length(gene) < 300)[1] == chr.genes[2]
        @test length(@genes(chr, get(gene, :locus_tag, "") == "")) == 3
    end

    @testset "Broadcast" begin
        # Broadcasted assignment on existing property
        chr.genes[3:4].gene .= "AXL2P"
        @test all(chr.genes[3:4].gene .== "AXL2P")
        # Broadcasted assignment on previously missing property
        chr.genes[3:4].newproperty .= true
        @test all(chr.genes[3:4].newproperty .== true)
        # Broadcasted assignment with @genes
        @genes(chr, :feature == "gene").newproperty .= false
        @test all(chr.genes[[3,5]].newproperty .== false)
    end

    @testset "Locus" begin
        locus = Locus(1:1, '.', true, true, UnitRange{Int}[])
        @test Locus() == locus
        @test chr.genes[2].locus < chr.genes[4].locus
        @test chr.genes[2].locus == chr.genes[2].locus
        @test iscomplement(chr.genes[2]) == false
        @test iscomplement(chr.genes[5]) == true
    end

    seq = dna"atgtccatatacaacggtatctccacctcaggtttagatctcaacaacggaaccattgccgacatgagacagttaggtatcgtcgagagttacaagctaaaacgagcagtagtcagctctgcatctgaagccgctgaagttctactaagggtggataacatcatccgtgcaagaccaagaaccgccaatagacaacatatgtaa"
    @test sequence(chr.genes[2]) == seq
    @test length(chr.genes[2]) == length(seq)

    @testset "Chromosome" begin
        chr = Chromosome()
        @test chr.name == ""
        @test chr.sequence == dna""
        @test chr.header == ""
        @test chr.genes == Gene[]
        @test names(chr.genedata) == [:feature, :locus]
    end
end
