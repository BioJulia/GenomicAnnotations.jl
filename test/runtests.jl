using GenomicAnnotations
using Base.Test

@testset "GenomicAnnotations" begin
    chr = readgbk("example.gbk")

    @testset "Gene properties" begin
        @test chr.genes[2].locus_tag == "tag01"
        @test (chr.genes[2].locus_tag = "tag01") == "tag01"
        @test begin
            chr.genes[2].test = "New column"
            chr.genes[2].test == "New column" && all(ismissing, chr.genes[[1,3,4,5,6]])
        end
    end

    @testset "Adding/removing genes" begin
        addgene!(chr, "CDS", 300:390, locus_tag = "tag04")
        @test chr.genes[end].locus_tag == "tag04"
        delete!(chr.genes[end])
        @test chr.genes[end].locus_tag == "tag03"
    end

    @testset "@genes" begin
        @test @genes(chr, :feature == "CDS") == chr.genes[[2,4,6]]
        @test @genes(chr, :complement) == chr.genes[[5,6]]
        @test @genes(chr, :feature == "CDS", ! :complement) == chr.genes[[2,4]]
        @test @genes(chr, length(gene) < 300) == chr.genes[2]
    end

    seq = "atgtccatatacaacggtatctccacctcaggtttagatctcaacaacggaaccattgccgacatgagacagttaggtatcgtcgagagttacaagctaaaacgagcagtagtcagctctgcatctgaagccgctgaagttctactaagggtggataacatcatccgtgcaagaccaagaaccgccaatagacaacatatgtaa"
    @test genesequence(chr.genes[2]) == seq
    @test length(chr.genes[2]) == length(seq)
end
