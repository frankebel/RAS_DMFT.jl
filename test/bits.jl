using RAS_DMFT
using Fermions
using Test

@testset "mask" begin
    @testset "mask_fe" begin
        f, e = mask_fe(UInt64, 1, 1, 1)
        @test f === UInt64(0b010_010)
        @test e === UInt64(0b100_100)
        f, e = mask_fe(UInt64, 4, 8, 9)
        @test f === UInt64(0b000000000111111110000_000000000111111110000)
        @test e === UInt64(0b111111111000000000000_111111111000000000000)
        f, e = mask_fe(BigMask{2}, 10, 11, 12)
        @test f === BigMask(
            (
                UInt64(0b0000000000111111111110000000000_000000000000111111111110000000000),
                UInt64(0b00),
            )
        )
        @test e === BigMask(
            (
                UInt64(0b1111111111000000000000000000000_111111111111000000000000000000000),
                UInt64(0b11),
            )
        )

        # wrong arguments
        @test_throws ArgumentError mask_fe(UInt64, 0, 0, 33)
        @test_throws ArgumentError mask_fe(UInt64, -1, 1, 1)
        @test_throws ArgumentError mask_fe(UInt64, 1, -1, 1)
        @test_throws ArgumentError mask_fe(UInt64, 1, 1, -1)

        @test_throws MethodError mask_fe(Int64, 1, 1, 1)
    end # mask_fe

    @testset "slater_start" begin
        @test slater_start(UInt64, 0b0000, 1, 2, 3, 4) ===
            UInt64(0b00000_111_00_1_00_0000_111_00_1_00)
        @test slater_start(BigMask{2}, 0b1001, 10, 11, 12, 13) ===
            BigMask(
            (
                0b0000_1111111111_10_0000000000000_111111111111_00000000000_1111111111_01,
                0b00000000000000000000000000000000_0000000000000_111111111111_0000000,
            )
        )

        # wrong arguments
        # 31 + 2 (bi) = 33 bits which is too much
        @test_throws ArgumentError slater_start(UInt64, zero(UInt8), 31, 0, 0, 0)
        @test_throws ArgumentError slater_start(UInt64, zero(UInt8), 0, 31, 0, 0)
        @test_throws ArgumentError slater_start(UInt64, zero(UInt8), 0, 0, 31, 0)
        @test_throws ArgumentError slater_start(UInt64, zero(UInt8), 0, 0, 0, 31)
        # negative values
        @test_throws ArgumentError slater_start(UInt64, zero(UInt8), -1, 0, 0, 0)
        @test_throws ArgumentError slater_start(UInt64, zero(UInt8), 0, -1, 0, 0)
        @test_throws ArgumentError slater_start(UInt64, zero(UInt8), 0, 0, -1, 0)
        @test_throws ArgumentError slater_start(UInt64, zero(UInt8), 0, 0, 0, -1)
    end # slater_start
end # mask
