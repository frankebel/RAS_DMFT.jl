using Aqua
using RAS_DMFT
using Test

@testset "Aqua" begin
    Aqua.test_all(RAS_DMFT)
end
