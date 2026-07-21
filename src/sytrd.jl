# sytrd: Reduce a real symmetric matrix to tridiagonal form.

using LinearAlgebra.BLAS: @blasfunc
using LinearAlgebra.LAPACK: chklapackerror
using LinearAlgebra:
    BlasInt,
    DimensionMismatch,
    checksquare,
    chkstride1,
    dot,
    libblastrampoline

for (sytrd, orgtr, elty) in ((:dsytrd_, :dorgtr_, :Float64), (:ssytrd_, :sorgtr_, :Float32))
    @eval begin
        function sytrd!(uplo::AbstractChar, A::AbstractMatrix{$elty})
            chkstride1(A)
            n = checksquare(A)
            d = similar(A, $elty, n)
            e = similar(A, $elty, n - 1)
            tau = similar(A, $elty, n - 1)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i in 1:2  # first call returns lwork as work[1]
                ccall(
                    (@blasfunc($sytrd), libblastrampoline),
                    Cvoid,
                    (
                        Ref{UInt8},
                        Ref{BlasInt},
                        Ptr{$elty},
                        Ref{BlasInt},
                        Ptr{$elty},
                        Ptr{$elty},
                        Ptr{$elty},
                        Ptr{$elty},
                        Ref{BlasInt},
                        Ptr{BlasInt},
                        Clong,
                    ),
                    uplo,
                    n,
                    A,
                    max(1, stride(A, 2)),
                    d,
                    e,
                    tau,
                    work,
                    lwork,
                    info,
                    1,
                )
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            return d, e, tau
        end

        function orgtr!(
                uplo::AbstractChar, A::AbstractMatrix{$elty}, tau::AbstractArray{$elty}
            )
            chkstride1(A)
            n = checksquare(A)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i in 1:2  # first call returns lwork as work[1]
                ccall(
                    (@blasfunc($orgtr), libblastrampoline),
                    Cvoid,
                    (
                        Ref{UInt8},
                        Ref{BlasInt},
                        Ptr{$elty},
                        Ref{BlasInt},
                        Ptr{$elty},
                        Ptr{$elty},
                        Ref{BlasInt},
                        Ptr{BlasInt},
                        Clong,
                    ),
                    uplo,
                    n,
                    A,
                    max(1, stride(A, 2)),
                    tau,
                    work,
                    lwork,
                    info,
                    1,
                )
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            return
        end
    end
end
