#################################################################################################################
#
#  This is a module used to approximate the values of the kernel of the fractional Laplacian on S^2 and other
#  subordinate processes of fractional order ν
#
#################################################################################################################
################################################################################################################
#
# We should start with the following command to compute the Fourier coefficients of the sampling functions
#
# P is the number of points in the sample.
# P = 10  
# a,m = pre_sample(P)
#
# Then, a and m can be used several times to sample with different values of ν and ω
# 
# Then, we can sample the values of a subordinate kernel with weight Cν t^(-2-ν)ω(t,ν) at P equispaced points from 0 to π/2
# with the command.
#
# sample_subordinate(νs,ω,a,m)
#
# Here, νs is a list of values of ν, and ω is a function that must depend on t and ν. The constant Cν
# is chosen so that ω=1 corresponds to the kernel of the fractional Laplacian.
#
# If we only want to sample the values of the kernel of the fractional Laplacian, instead of calling sample_subordinate
# with ω=1, it is more efficient to call
#
# sample_fractional_laplacian(νs,a,m)
#
##############################################################################################################




module Subordinate

export evenY, Cnu
export eigenvals_from_kernel, subordinate_eigenvals
export pre_sample, sample_subordinate, sample_fractional_laplacian

using SphericalHarmonics
using QuadGK
using SpecialFunctions

# We precompute all the coefficients of the spherical harmonics
sphericalcoefficients = SphericalHarmonics.compute_coefficients(10000, 0)

function Yl0(θ::Real, index::Integer)
    return SphericalHarmonics.sphericalharmonic(θ,0.,l=index,m=0,SHType=SphericalHarmonics.RealHarmonics(),coeff=sphericalcoefficients)
end

function evenY(θ::Real, index::Integer)
    return Yl0(θ,2*index)
end

function Cnu(ν::Real)
    return 2^ν * gamma(ν/2+1) / π / abs(gamma(-ν/2)) / 2
end

######################################################################################
#################### Generic kernels #################################################
######################################################################################
######################################################################################


# This function starts from a kernel b(θ) and returns the eigenvalues of the
# corresponding integro-differential operator for Y_1, Y_2, Y_3, ..., Y_N

function eigenvals_from_kernel(b::Function, N::Integer)
    lambda = zeros(N)
    for index in 1:N
        f(θ) = evenY(θ,index)
        g(θ) = 2π * (f(0)-f(θ)) * b(θ) * sin(θ)
        Bf = quadgk(g,0,π)[1]
        lambda[index] = Bf / f(0)
    end
    return lambda
end

# The following function computes the eigenvalues of a subordinate kernel with weight ω
# it returns a list of eigenvalues corresponding to Y_1, Y_2, ..., Y_N
function subordinate_eigenvals(ω::Function, N::Integer)
    lambda = zeros(N)
    verymuch = 20000.
    for index in 1:N
        λ = 2index*(2index+1)
        g(t) = 1-exp(-t)
        integrand(t) = g(t)*ω(t/λ)/λ
        lambda[index] = quadgk(integrand,0,verymuch)[1]
    end
    return lambda
end

# The following function computes the eigenvalues of a subordinate kernel with weight ω(t) * t^{-1-ν/2} / factor
# assuming that ω(t) = 1 for t large. It is more accurate than the previous implementation under that assumption
# it returns a list of eigenvalues corresponding to Y_1, Y_2, ..., Y_N
# The should be the eigenvalues of (-Δ)^{ν/2} when ω=1.
function subordinate_eigenvals(ω::Function, ν::Real, N::Integer)
    lambda = zeros(N)
    verymuch = 100.
    g(t) = (1-exp(-t))
    ω1(t,λ) = ω(t/λ)*t^(-1-ν/2)*λ^(ν/2)

    for index in 1:N
        λ = 2index*(2index+1)
        #ω1(t) = ω(t/λ)*t^(-1-ν/2)*λ^(ν/2)
        factor = subordinate_factor(ν)
        #integrand(t) = g(t)*ω1(t)
        term1 = quadgk(t -> g(t)*ω1(t,λ),0,verymuch)[1]
        
        ###### t = ((1-ν/2)*s)^(1/(1-ν/2))
        #s1 = (10^(1-ν/2))/(1-ν/2)
        #int1(s) = ω( ((1-ν/2)*s)^(1/(1-ν/2)) /λ) * (g( ((1-ν/2)*s)^(1/(1-ν/2)) ) / ((1-ν/2)*s)^(1/(1-ν/2)))
        #term11 = λ^(ν/2)*quadgk(int1,0.,s1)[1]
        #term12 = quadgk(integrand,10,verymuch)[1]
        #term1 = term11 + term12

        term2 = verymuch^(-ν/2)*λ^(ν/2)*2/ν*ω(verymuch/λ)
        lambda[index] = (term1+term2)/factor
    end
    return lambda
end


# The following function computes the factor so that factor*t^{-1-ν/2} is the subordinate
# function for the fractional Laplacian (-Δ)^{ν/2}
function subordinate_factor(ν::Real)
    verymuch = 100.
    return quadgk(s->(1-exp(-s))*s^(-1-ν/2),0,verymuch)[1] + 2/ν * verymuch^(-ν/2)
end

# The following functions returns the coefficients a so that we reproduce the integro-differential kernel 
# on S^2 from the subordinate weight t^(-1-2ν) ω(t) / subordinate_factor
# It should return the coefficients of the fractional Laplacian when ω=1, but the numerical error is quite significant
# The coefficients correspond to Y_0, Y_1, ..., Y_{N-1}

function subordinate_coefficients(ν::Real, ω::Function, N::Integer)
    factor = subordinate_factor(ν)
    ω1(t) = ω(t)
    lambda = subordinate_eigenvals(ω1,ν,N)
    return compute_b_coefficients(ν,lambda)
end

# This one simply evaluates a function as a series of even spherical harmonics

function harmonic_series(θ::Real, a::Vector{<:Real})
    N = length(a)
    res = 0.
    for index in N:-1:1
        res += a[index]*evenY(θ,index-1)
    end
    return res
end

################################# Sampling method ########################################

# We implement our own integration rule, to keep things under control.

function simpsons(f::Function,a,b;N=10000)
    res = (f(a)+f(b))/6
    for i in 1:2*N-1
        x = a+i*(b-a)/(2N)
        if (mod(i,2)==1)
            res += 2/3 * f(x)
        else
            res += 1/3 * f(x)
        end
    end
    return res * (b-a) / N
end

# Returns the Fourier coefficients from a given function on S^2
function series_coefficients(f::Function, N::Integer; precision=10000, support=(0,π))
    a = zeros(N)
    for ℓ in 1:N
        a[ℓ] = simpsons(θ -> 2π*sin(θ)*evenY(θ,ℓ-1)*f(θ),support[1],support[2];N=precision)
    end
    return a
end

function fractional_laplacian(f::Function, ν::Real, N::Integer; tol=0.)
    if tol==0.
        a = series_coefficients(f,N)
    else
        a = series_coefficients(f,N,tol)
    end
    return fractional_laplacian(a,ν)
end

function fractional_laplacian(a::Vector{<:Real}, ν::Real)
    N = length(a)
    b = similar(a)
    for ℓ in 1:N
        b[ℓ] = a[ℓ]*(2(ℓ-1)*(2ℓ-1))^(ν/2)
    end    
    return harmonic_series(0.,b)
end

function fourier_multiplier(a::Vector{<:Real}, m::Vector{<:Real})
    @assert length(a) == length(m)
    N = length(a)
    b = similar(a)
    for ℓ in 1:N
        b[ℓ] = a[ℓ]*m[ℓ]
    end    
    return harmonic_series(0.,b)
end

function pre_sample(points::Integer; N::Integer=0, gg=0)
    P = points
    mass = zeros(P)
    if (gg==0)
        gg = P^2
    end
    prec = div(gg^2 , 20)
    if N==0
        N = 30*points
    end
    a = zeros(P,N)

    Threads.@threads for i in 1:P
        θ0 = π/2 * (P-i) / (P+1)

        # For the first couple of points, we increase the accuracy.
        gg2 = gg
        prec2 = prec
        if (i==1)
            gg2 = 4*gg
            prec2 = 10*prec
        end

        if i in [2,3,4]
            gg2 = 2*gg
            prec2 = 4*prec
        end


        if (gg2 < 10.)
        	# This is a legacy condition that should only be executed if the computation is intended to be very fast
        	# and inaccurate.
            f1(θ) = exp(-gg2*(cos(2θ-2θ0)+1))+exp(-gg2*(cos(2θ+2θ0)+1))
            mass[i] = simpsons(θ -> 4π*f1(θ)*sin(θ),0,π/2)
            a[i,:] = series_coefficients(f1,N,precision=prec2)
        else
        	# This is a function centered at π/2 - θ0 so that
        	# \int f(θ) (θ-π/2+θ0) dθ = 0 
        	# \int f(θ) (θ-π/2+θ0)^2 dθ = 0
        	# Making the second moment equal to zero should decrease the final error term by an order of magnitude. 
            f2(θ) = 2sqrt(2)*exp(-4*gg2*(π/2 - θ0 - θ)^2) - exp(-2*gg2*(π/2 - θ0 - θ)^2)
            supp = sqrt(20/gg2)
            @assert π/2 - θ0 - supp > 0
            @assert π/2 - θ0 + supp < π
            mass[i] = simpsons(θ -> 4π*f2(θ)*sin(θ),π/2 - θ0-supp,π/2 - θ0+supp)
            prec2 = 400
            a[i,:] = 2*series_coefficients(f2,N,precision=prec2,support=(π/2 - θ0-supp,π/2 - θ0+supp))
        end
    end
    return a, mass
end    

function sample_subordinate(ν::Vector{<:Real}, ω::Function, a::Matrix{<:Real}, mass::Vector{<:Real})
    P = length(mass)
    N = size(a)[2]
    M = length(ν)    
    res = zeros(M,P)

    Threads.@threads for n in 1:M
        m1 = subordinate_eigenvals(t->ω(t,ν[n]), ν[n], N-1)
        m = vcat([0],m1)
    
        for i in 1:P
            res[n,i] = -fourier_multiplier(a[i,:],m)/mass[i]
        end
    end
    return res
end

function sample_fractional_laplacian(ν::Vector{<:Real}, a::Matrix{<:Real}, mass::Vector{<:Real})
    P = length(mass)
    N = size(a)[2]
    M = length(ν) 
    res = zeros(M,P)

    Threads.@threads for n in 1:M
        for i in 1:P
            res[n,i] = -fractional_laplacian(a[i,:],ν[n])/mass[i]
        end
    end
    return res
end

end