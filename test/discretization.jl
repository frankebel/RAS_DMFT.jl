using RAS_DMFT
using Test

@testset "discretization" begin
    @testset "discretize_to_grid" begin
        W = [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        f = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0] # not same step size
        # 1 pole
        grid = [0.0]
        foo = discretize_to_grid(f, W, grid)
        @test locations(foo) == grid
        @test locations(foo) !== grid
        @test weights(foo) == [15.125 / π]
        # 1 pole < W
        foo = discretize_to_grid(f, W, [-10.0])
        @test locations(foo) == [-10.0]
        @test weights(foo) == [15.125 / π]
        # 1 pole > W
        foo = discretize_to_grid(f, W, [10.0])
        @test locations(foo) == [10.0]
        @test weights(foo) == [15.125 / π]
        # 2 asymmetric poles
        foo = discretize_to_grid(f, W, [-3.0, 5.0])
        @test locations(foo) == [-3.0, 5.0]
        @test weights(foo) == [6.0 / π, 9.125 / π]
    end # discretize_to_grid
end # discretization
