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
function uofr(r::Real,p::Real,q::Real,ψ0::Real=1/4)
    return 1 - p^2 / r^2 - 4ψ0 / r^(q-1)
end

function ∂u∂r(r::Real,p::Real,q::Real,ψ0::Real=1/4)
    return 2*p^2/r^3 + 4*(q-1)*ψ0 / r^q
end

# This is the function that returns r for any given u
# We use it to compute r0
# We invert the function uofr using Newton's method
function rofu(u::Real,p::Real,q::Real,ψ0::Real=1/4)
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

function integrand3(u::Real, p::Real, q::Real, r0::Real, ψ0::Real=1/4)
    rr1 = 1-u^2
    #rr2 = r0^2/p^2 - u^2 - u^(q-1) / (r0^(q-3) * p^2)
    rr2 = 1 - u^2 + (1 - u^(q-1)) / (r0^(q-3) * p^2)
    if ((rr1<=0) || (rr2<=0)) return 0. end
    term1 = 2 / sqrt(rr1)
    term2 = 2 / sqrt(rr2)
    return term1 - term2
end

function θofp(p::Real, q::Real, ψ0::Real=1/4)
    r0 = rofu(0,p,q,ψ0)
    f(u) = integrand3(u,p,q,r0)
    term1 = quadgk(f,0,1)[1]
    return  term1
end

# We differentiate \theta(p) numerically with an incremental quotient
function ∂θ∂p(p::Real, q::Real, ψ0::Real=1/4)
    h = 0.001
    pos = θofp(p+h,q,ψ0) + θofp(p+2h,q,ψ0)
    neg = θofp(p-h,q,ψ0) + θofp(p-2h,q,ψ0)
    return (pos-neg)/(6*h)
end

function bisection(f::Function,x0::Real,x1::Real,tolerance::Real=1e-15)
    counter = 0
    while abs(x1-x0) > tolerance
        x2 = (x1+x0)/2
        if f(x2)>=0 
            x1 = x2
        else
            x0 = x2
        end
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
function pofθ(θ::Real, q::Real, ψ0::Real=1/4)
    ff = r -> θ - θofp(r,q,ψ0)
    p0 = 0.
    # I am not 100% sure this is always large enough.
    # Choosing a larger value for p1 makes this method much slower
    p1 = 2*(1+q)*ψ0/θ^(1/(q-1)) 
    tol = 1.e-10
    return bisection(ff,p0,p1,tol)
end

# We are now ready to write the kernel.
# The factor1 is the only difference with the computation on S^1

function kernelb(θ::Real, q::Real, ψ0::Real=1/4; dimension=3)
    p = pofθ(θ,q,ψ0)
    factor1 = p / sin(θ)
    der = -∂θ∂p(p,q,ψ0)
    return (factor1)^(dimension-2)/der
end

# Since we only care about even functions, we may symmetrize the kernel b
function symb(θ::Real, q::Real, ψ0::Real=1/4; dim=3)
    return kernelb(θ,q,ψ0,dimension=dim) + kernelb(π-θ,q,ψ0,dimension=dim)
end

function νofq(q::Real; dimension=3)
	return (dimension-1)/(q-1)
end

function qofν(ν::Real; dimension=3)
	return 1 + (dimension-1)/ν
end

function Cq(q::Real; dimension=3)
   return (sqrt(π)*Γ(q/2) / Γ(q/2-1/2))^((dimension-1)/(q-1))/(q-1)
end


end