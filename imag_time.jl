using SparseArrays
using LinearAlgebra
using PyPlot

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

    E = sum(1/(2*h).*(ϕ[2:end]-ϕ[1:end-1]).^2) +
        sum(h.*V.(x).*(ϕ.^2) + h*β/2 .*ϕ.^4) +
        1/(2*h)*(ϕ[1]^2 + ϕ[end]^2)/(2*h)
    μ = E + sum(h*β/2 .*ϕ.^4)
    x̄ = zeros(M+1)
    ϕ₁= zeros(M+1)
    x̄[1] = a
    x̄[2:end-1] = x
    x̄[end] = b
    ϕ₁[2:end-1] = ϕ
    return x̄, ϕ₁, E, μ
end

function imag_time_2d(β,a,b,M,k,n,V,ϕ₀)
    h = (b-a)/M
    x̄ = LinRange(a,b,M+1)
    x = x̄[2:end-1]
    Δ₀ = sparse(laplacian_1d(M-1,h))
    I₀ = sparse(I,M-1,M-1)
    Δ = kron(Δ₀,I₀) + kron(I₀,Δ₀)
    V₀ = V.(x,x')
    V̄ = sparse(diagm(vec(V₀)))
    ϕ = vec(ϕ₀.(x,x'))
    Id = sparse(I,(M-1)^2,(M-1)^2)

    for i = 1:n
        A = V̄ + sparse(diagm(β.*(ϕ).^2))
        ϕ̃ = -((k/2)*Δ - k*A - Id)\ϕ
        ϕ = ϕ̃ ./ (h*norm(ϕ̃))
    end

    Φ = zeros(M+1,M+1)
    Φ[2:end-1,2:end-1] = reshape(ϕ,(M-1,M-1))

    if sum(Φ) < 0
        Φ = -Φ
    end

    E = 1/2 * sum((Φ[2:end,:]-Φ[1:end-1,:]).^2) +
        1/2 * sum((Φ[:,2:end]-Φ[:,1:end-1]).^2) +
        h^2 * sum(V₀.*((Φ[2:end-1,2:end-1]).^2)) +
        β*h^2/2 * sum(Φ.^4)
    μ = E + β*h^2/2 * sum(Φ.^4)

    return x̄, Φ, E, μ
end

function TF_approx(V,β,x)
    μ = 1/2*(3*β/2)^(2/3)
    if V(x) > μ
        return 0
    end
    return sqrt((μ-V(x))/β)
end

function compare(β,a,b,M,V)
    x, ϕ₁, E₁, μ = imag_time_1d(β,a,b,M,0.1,200,V,x->(x-a)*(x-b))
    ϕ₂ = TF_approx.(V,β,x)
    E₂ = 3/10*(3*β/2)^(2/3)
    clf()
    plot(x, ϕ₁, label = "imaginary time")
    plot(x, ϕ₂, label =  "Thomas-Fermi")
    legend(loc="best")
    savefig("β=$(β).png")
    return E₁,E₂
end

function example()
    V(x) = x^2/2 + 12*sin(x/2)^2
    ϕ(x) = 16^2 - x^2

    clf()
    for β in  [0,5,50,100,500,1000,2000]
        x, ϕ₁, E, μ = imag_time_1d(β,-16,16,128,0.1,200,V,ϕ)
        E = round(E,digits=3)
        μ = round(μ,digits=3)
        plot(x, ϕ₁, label = "β=$(β), μ=$(μ), E=$(E)")
    end
    legend(loc=2)
    savefig("imag_time.png")
end
