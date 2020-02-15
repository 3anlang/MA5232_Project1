using FFTW
using SparseArrays
using LinearAlgebra
using Plots

function laplacian_1d(n,h)
    return 1.0/h^2 * Tridiagonal(
    fill( 1.0,n-1),
    fill(-2.0,n),
    fill( 1.0,n-1))
end

function imag_time_1d(β,a,b,M,k,n,V,ϕ₀)
    h = (b-a)/M
    x̄ = LinRange(a,b,M+1)
    x = x̄[2:end-1]
    Δ = sparse(laplacian_1d(M-1,h))
    Id = sparse(I,M-1,M-1)
    V̄ = sparse(diagm(V.(x)))
    ϕ = ϕ₀.(x)

    for i = 1:n
        A = V̄ + sparse(diagm(β.*(ϕ).^2))
        ϕ̃ = -((k/2)*Δ - k*A - Id)\ϕ
        ϕ = ϕ̃ ./ (sqrt(h)*norm(ϕ̃))
    end

    if sum(ϕ) < 0
        ϕ = -ϕ
    end

    ϕ₁= zeros(M+1)
    ϕ₁[2:end-1] = ϕ
    return x̄, ϕ₁, h
end

function tssp_1d(β,a,b,M,k,n,v,ϕ₀,α)
    h = (b-a)/(M-1)
    x = LinRange(a,b,M)
    ϕ = ϕ₀[1:M]
    V = v.(x)
    Φ = zeros(Complex{Float64},n,M)
    μ = 2*im*π/(M*h).*(0:M-1)
    for i = 1:n
        Φ[i,:] = ϕ
        ϕ¹ = exp.((-im*k/2)*(V + β.*(abs.(ϕ).^2))).*ϕ
        Ψ = fft(ϕ¹)
        Ψ¹ = exp.(im*k/2 .*(μ.^2)).*Ψ
        ϕ = ifft(Ψ¹)
        ϕ =  exp.((-im*k/2)*(V + β.*(abs.(ϕ).^2))).*ϕ
    end

    max = maximum(abs.(Φ))
    anim = @animate for i=1:n
        plot(x, abs.(Φ[i,:]).^2, label = "t = $(round(i*k,digits=4))")
        ylims!((0.0,(max^2)*5/4))
    end
    gif(anim, "1D α=$(α) β=$(β).gif", fps = 5)
    return Φ,h,x
end

function dyna_1d(β,a,b,M,k,n,V,α,x₀)
    f = x -> (x-a)*(b-x)
    x, ϕ₀, h = imag_time_1d(β,a+x₀,b+x₀,M,k,n,V,f)
    x = LinRange(a,b,M+1)
    ψ₀ = ϕ₀.*exp.(im*α*x)

    Φ,h,x = tssp_1d(β,a,b,M,k,n,V,ψ₀,α)
    x̄ = zeros(n)
    σ = zeros(n)
    for i = 1:n
        x̄[i] = h*sum(x.*abs.(Φ[i,:]).^2)
        σ[i] = sqrt(h*sum(x.^2 .*abs.(Φ[i,:]).^2))
    end
    plot(1:n,x̄,label = "x(t)")
    savefig("Q9(1) α=$(α) β=$(β).png")
    plot(1:n,σ,label = "sigma(t)")
    savefig("Q9(2) a=$(α) β=$(β).png")
    plot(1:n,σ.^2-x̄.^2,label = "sigma^2(t)-x^2")
    savefig("Q9(3) a=$(α) β=$(β).png")
    return x̄,σ
end
