using Test
using Random
using GeneRegulatorySystems
using GeneRegulatorySystems.Models: Models, FlatState

@testset "reproducibility" begin
    @testset "branched V1 schedule is parent-seed-sensitive and reproducible" begin
        spec = """
        {
            "step": {
                "do": {"{regulation/v1}": {"genes": [{"\$": ["defaults", "gene"]}]}},
                "step": [
                    {"{add}": {"\$": ["defaults", "bootstrap"], "1.proteins": 100}},
                    {"branch": true, "each": {"length": 3}, "as": "i",
                     "step": [{"to": 100.0}]}
                ]
            }
        }
        """
        f! = Models.parse(spec; seed = "")
        run_with(s) = let x = FlatState(randomness = Xoshiro(s))
            [sort(collect(FlatState(b).counts)) for b in f!(x).branches]
        end

        a = run_with(UInt64(42))
        b = run_with(UInt64(42))
        c = run_with(UInt64(43))

        @test a == b                  # reproducible
        @test length(unique(a)) == 3  # branches independent
        @test a != c                  # parent-seed sensitive
    end
end
