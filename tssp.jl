using FFTW
using SparseArrays
using LinearAlgebra
using Plots

function tssp_1d(β,a,b,M,k,n,v,ϕ₀)
    h = (b-a)/(M-1)
    x = LinRange(a,b,M)
    ϕ = ϕ₀[1:M]
    V = v.(x)
    Φ = zeros(Complex{Float64},n,M)
    μ = 2*im*π/(M*h).*(0:M-1)
    for i = 1:n
        Φ[i,:] = ϕ
        ϕ¹ = exp.((-im*k/2)*(V + β.*abs.(ϕ).^2)).*ϕ
        Ψ = fft(ϕ¹)
        Ψ¹ = exp.(im*k/2 .*μ.^2).*Ψ
        ϕ = ifft(Ψ¹)
        ϕ =  exp.((-im*k/2)*(V + β.*abs.(ϕ).^2)).*ϕ
    end

    anim = @animate for i=1:n
        plot(x, abs.(Φ[i,:]).^2, label = "t = $(round(i*k,digits=4))")
        ylims!((0.0,1.0))
        title!("1D Bose-Einstein Condesation Dynamics")
    end
    gif(anim, "1D.gif", fps = 5)
    return Φ,h,x
end

function tssp_2d(β,a,b,M,k,n,v,ϕ₀)
    h = (b-a)/(M-1)
    x = LinRange(a,b,M)
    ϕ = ϕ₀[1:M,1:M]
    V = v.(x,x')
    Φ = zeros(Complex{Float64},n,M,M)
    μ = 2*im*π/(M*h).*(0:M-1)
    for i = 1:n
        Φ[i,:,:] = ϕ
        ϕ¹ = exp.((-im*k/2)*(V + β.*abs.(ϕ).^2)).*ϕ
        Ψ = fft(ϕ¹)
        Ψ¹ = exp.(im*k/2 .*kron(μ.^2,(μ.^2)')).*Ψ
        ϕ = ifft(Ψ¹)
        ϕ = exp.((-im*k/2)*(V + β.*abs.(ϕ).^2)).*ϕ
    end

    function f(x,y,Φ)
        i = convert(Int64, round((x-a)*(M-1)/(b-a)+1,digits=0))
        j = convert(Int64, round((y-a)*(M-1)/(b-a)+1,digits=0))
        return abs(Φ[i,j])^2
    end

    anim = @animate for i=1:n
        plot(x,x, (x,y)->f(x,y,Φ[i,:,:]),st=wireframe)
        zlims!((0,0.12))
        title!("2D Bose-Einstein Condesation Dynamics")
    end
    gif(anim, "2D.gif", fps = 5)
    return Φ,h,x
end
