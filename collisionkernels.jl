##############################################################################################
#
# kernelb(θ,q) 
# returns the collision kernel corresponding to power-law potentials with exponent q
# This is the textbook value of the kernel. We should have q>2.
#
# Since we only care about even functions on the sphere, the convenient function symb is the
# symmetrization of the kernel
# symb(θ,q) = kernelb(θ,q) + kernelb(π-θ,q)
#
# The order of the integro-differential operator ν is obtained with 
# νofq(q)
# and inversely with
# qofν(ν)
#
# The assymptotic behavior of the kernel near θ=0 should be ≈Cq θ^(-2-ν)
# The value of Cq is obtained with 
# Cq(q)
#
# This file has been tested on Julia 1.10.4
#
##############################################################################################

module CollisionKernels
export kernelb, symb, νofq, qofν, Cq

using SpecialFunctions
using QuadGK
Γ = gamma

#This is the function that returns u for any given value of r
function uofr(r::Real,p::Real,q::Real,ψ0::Real=1)
    return 1 - p^2 / r^2 - 4ψ0 / r^(q-1)
end

function ∂u∂r(r::Real,p::Real,q::Real,ψ0::Real=1)
    return 2*p^2/r^3 + 4*(q-1)*ψ0 / r^q
end

# This is the function that returns r for any given u
# We use it to compute r0
# We invert the function uofr using Newton's method
function rofu(u::Real,p::Real,q::Real,ψ0::Real=1)
    @assert 0. <= u <= 1.
    r1 = 0.
    r2 = max(p,(4ψ0)^(1/q))
    tol = 1e-7
    while abs(r1-r2) > tol
        r3 = r2 + (u-uofr(r2,p,q,ψ0))/ ∂u∂r(r2,p,q,ψ0)
        r1 = r2
        r2 = r3
        # println(r2," : ",uofr(r2,p,q,ψ0))
    end
    return r2
end

#function integrand(r::Real, p::Real, q::Real, ψ0::Real=1)
#    den = 1. - p^2 / r^2 - 4ψ0 / r^(q-1)
#    tol = 1.e-7
#    if den <= tol return 0. end
#    return 1/r^2 / sqrt(den)
#end

function integrand2(r::Real, p::Real, q::Real, ψ0::Real=1)
    α = 4ψ0 / (p^2 * r^(q-3))
    den = 1. - p^2 / r^2 - 4ψ0 / r^(q-1)
    tol = 1.e-7
    if den <= tol return 0. end
    num = 1 + (q-1)*α/2 - sqrt(1+α)
    return 2p / r^2 / sqrt(den) * num / sqrt(1+α)
end

# This is the function that computes θ for any given p
# function θofp(p::Real, q::Real, ψ0::Real=1)
#     f(r) = integrand(r,p,q,ψ0)
#     r0 = rofu(0,p,q,ψ0)
#     return π - 2*p*quadgk(f,r0,Inf)[1]
# end

function θofp(p::Real, q::Real, ψ0::Real=1)
    f(r) = integrand2(r,p,q,ψ0)
    r0 = rofu(0,p,q,ψ0)
    # verymuch = 10000*(100+p)
    term1 = quadgk(f,r0,Inf)[1]  # For some reason, the integral takes a very long time when p>>1
    #term2 = 4ψ0*(q-2)*verymuch^(-q+2)/(p*(q-2))
    return  term1 #+ term2
end

# We differentiate \theta(p) numerically with an incremental quotient
function ∂θ∂p(p::Real, q::Real, ψ0::Real=1)
    h = 0.001
    pos = θofp(p+h,q,ψ0) + θofp(p+2h,q,ψ0)
    neg = θofp(p-h,q,ψ0) + θofp(p-2h,q,ψ0)
    return (pos-neg)/(6*h)
end

function bisection(f::Function,x0::Real,x1::Real,tolerance::Real=1e-15)
    counter = 0
    while abs(x1-x0) > tolerance
        # Replace the following line with the right formula for x2.
        # It should be the same one that you used for the secant method.
        x2 = (x1+x0)/2
        if f(x2)>=0 
            x1 = x2
        else
            x0 = x2
        end
        # println("RF: ",x0," , ",x1)
        counter = counter + 1
        if counter > 100000
            print("!")
            return (x0+x1)/2
        end
    end
    # println("We performed ", counter, " iterations.")
    return (x0+x1)/2
end

# This is the function that computes p for any given θ
function pofθ(θ::Real, q::Real, ψ0::Real=1)
    ff = r -> θ - θofp(r,q,ψ0)
    p0 = 0.
    # I am not 100% sure this is always large enough.
    # Choosing a larger value for p1 makes this method much slower
    p1 = 2*q*ψ0/θ^(1/(q-1)) 
    tol = 1.e-7
    return bisection(ff,p0,p1,tol)
end

# We are now ready to write the kernel.
# The factor1 is the only difference with the computation on S^1

function kernelb(θ::Real, q::Real, ψ0::Real=1.; dimension=3)
    p = pofθ(θ,q,ψ0)
    factor1 = p / sin(θ)
    der = -∂θ∂p(p,q,ψ0)
    return (factor1)^(dimension-2)/der
end

# Since we only care about even functions, we may symmetrize the kernel b
function symb(θ::Real, q::Real, ψ0::Real=1.)
    return kernelb(θ,q,ψ0) + kernelb(π-θ,q,ψ0)
end

function νofq(q::Real; dimension=3)
	return (dimension-1)/(q-1)
end

function qofν(ν::Real; dimension=3)
	return 1 + (dimension-1)/ν
end

function Cq(q::Real)
   ν = νofq(q)
   return (2*(q-2)*sqrt(π)*Γ(q/2-1) / Γ(q/2-1/2))^(2/(q-1))*ν/2
end

end