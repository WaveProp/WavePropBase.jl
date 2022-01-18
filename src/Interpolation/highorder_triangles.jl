# P3 for ReferenceTriangle
# source: https://www.math.uci.edu/~chenlong/iFEM/doc/html/dofP3doc.html
@fastmath function (el::LagrangeTriangle{10})(u)
    @assert u ∈ domain(el)
    λ₁ = 1-u[1]-u[2]
    λ₂ = u[1]
    λ₃ = u[2]
    ϕ₁ = 0.5*(3λ₁-1)*(3λ₁-2)*λ₁
    ϕ₂ = 0.5*(3λ₂-1)*(3λ₂-2)*λ₂
    ϕ₃ = 0.5*(3λ₃-1)*(3λ₃-2)*λ₃
    ϕ₄ = 4.5*λ₁*λ₂*(3λ₁-1) 
    ϕ₅ = 4.5*λ₁*λ₂*(3λ₂-1)
    ϕ₆ = 4.5*λ₃*λ₂*(3λ₂-1)
    ϕ₇ = 4.5*λ₃*λ₂*(3λ₃-1)
    ϕ₈ = 4.5*λ₁*λ₃*(3λ₃-1)
    ϕ₉ = 4.5*λ₁*λ₃*(3λ₁-1)
    ϕ₁₀ = 27*λ₁*λ₂*λ₃
    v = vals(el)
    return v[1]*ϕ₁ + v[2]*ϕ₂ + v[3]*ϕ₃ +
           v[4]*ϕ₄ + v[5]*ϕ₅ + v[6]*ϕ₆ +
           v[7]*ϕ₇ + v[8]*ϕ₈ + v[9]*ϕ₉ +
           v[10]*ϕ₁₀
end

@fastmath function jacobian(el::LagrangeTriangle{10,T}, u) where T
    @assert u ∈ domain(el)
    λ₁ = 1-u[1]-u[2]
    λ₂ = u[1]
    λ₃ = u[2]
    ∇λ₁ = SMatrix{1,2,eltype(T),2}(-1.,-1.)
    ∇λ₂ = SMatrix{1,2,eltype(T),2}(1.,0.)
    ∇λ₃ = SMatrix{1,2,eltype(T),2}(0.,1.)
    ∇ϕ₁ = (13.5*λ₁*λ₁-9λ₁+1)*∇λ₁
    ∇ϕ₂ = (13.5*λ₂*λ₂-9λ₂+1)*∇λ₂
    ∇ϕ₃ = (13.5*λ₃*λ₃-9λ₃+1)*∇λ₃
    ∇ϕ₄ = 4.5*((3*λ₁*λ₁-λ₁)*∇λ₂+λ₂*(6λ₁-1)*∇λ₁)
    ∇ϕ₅ = 4.5*((3*λ₂*λ₂-λ₂)*∇λ₁+λ₁*(6λ₂-1)*∇λ₂)
    ∇ϕ₆ = 4.5*((3*λ₂*λ₂-λ₂)*∇λ₃+λ₃*(6λ₂-1)*∇λ₂)
    ∇ϕ₇ = 4.5*((3*λ₃*λ₃-λ₃)*∇λ₂+λ₂*(6λ₃-1)*∇λ₃)
    ∇ϕ₈ = 4.5*((3*λ₃*λ₃-λ₃)*∇λ₁+λ₁*(6λ₃-1)*∇λ₃)
    ∇ϕ₉ = 4.5*((3*λ₁*λ₁-λ₁)*∇λ₃+λ₃*(6λ₁-1)*∇λ₁)
    ∇ϕ₁₀ = 27*(λ₁*λ₂*∇λ₃+λ₁*λ₃*∇λ₂+λ₃*λ₂*∇λ₁)
    v = vals(el)
    return v[1]*∇ϕ₁ + v[2]*∇ϕ₂ + v[3]*∇ϕ₃ +
           v[4]*∇ϕ₄ + v[5]*∇ϕ₅ + v[6]*∇ϕ₆ +
           v[7]*∇ϕ₇ + v[8]*∇ϕ₈ + v[9]*∇ϕ₉ +
           v[10]*∇ϕ₁₀
end

# P4 for ReferenceTriangle
# source: Silvester PP, Ferrari RL, Finite elements for electrical engineers (1990).
@fastmath function (el::LagrangeTriangle{15})(u)
    @assert u ∈ domain(el)
    x = u[1]
    y = u[2]
    ϕ₁ = 1/3*(-1+x+y)*(-1+2x+2y)*(-3+4x+4y)*(-1+4x+4y)
    ϕ₂ = 1/3*x*(-1+2x)*(-3+4x)*(-1+4x)
    ϕ₃ = 1/3*y*(-1+2y)*(-3+4y)*(-1+4y)
    ϕ₄ = -16/3*x*(-1+x+y)*(-1+2x+2y)*(-3+4x+4y)
    ϕ₅ = 4x*(-1+4x)*(-1+x+y)*(-3+4x+4y)    
    ϕ₆ = -16/3*x*(1-6x+8x^2)*(-1+x+y)
    ϕ₇ = 8/3*x*(-2+4x)*(-1+4x)*y
    ϕ₈ = 4x*(-1+4x)*y*(-1+4y)
    ϕ₉ = 8/3*x*y*(-2+4y)*(-1+4y)
    ϕ₁₀ = -16/3*y*(-1+x+y)*(1-6y+8y^2)
    ϕ₁₁ = 4y*(-1+x+y)*(-1+4y)*(-3+4x+4y)
    ϕ₁₂ = -16/3*y*(-1+x+y)*(-1+2x+2y)*(-3+4x+4y)
    ϕ₁₃ = 32x*y*(-1+x+y)*(-3+4x+4y)
    ϕ₁₄ = -32x*(-1+4x)*y*(-1+x+y)
    ϕ₁₅ = -32x*y*(-1+x+y)*(-1+4y)
    v = vals(el)
    return v[1]*ϕ₁ + v[2]*ϕ₂ + v[3]*ϕ₃ +
           v[4]*ϕ₄ + v[5]*ϕ₅ + v[6]*ϕ₆ +
           v[7]*ϕ₇ + v[8]*ϕ₈ + v[9]*ϕ₉ +
           v[10]*ϕ₁₀ + v[11]*ϕ₁₁ + v[12]*ϕ₁₂ +
           v[13]*ϕ₁₃ + v[14]*ϕ₁₄ + v[15]*ϕ₁₅
end

@fastmath function jacobian(el::LagrangeTriangle{15,T}, u) where T
    @assert u ∈ domain(el)
    x = u[1]
    y = u[2]
    L = SMatrix{1,2,eltype(T),2}
    ∇ϕ₁ = L(1/3*(-5+8x+8y)*(5+16x^2+4y*(-5+4y)+4x*(-5+8y)), 1/3*(-5+8x+8y)*(5+16x^2+4y*(-5+4y)+4x*(-5+8y)))
    ∇ϕ₂ = L(1/3*(-3+8x)*(1+4x*(-3+4x)), 0)  
    ∇ϕ₃ = L(0, 1/3*(-3+8y)*(1+4y*(-3+4y)))
    ∇ϕ₄ = L(-16/3*(-3+2x*(13+x*(-27+16x))+13y+72*(-1+x)*x*y+6*(-3+8x)*y^2+8y^3), -16/3*x*(13+24x^2+12y*(-3+2y)+12x*(-3+4y)))
    ∇ϕ₅ = L(4*(-1+2x+y)*(3-4y+32x*(-1+x+y)), 4x*(-1+4x)*(-7+8x+8y))  
    ∇ϕ₆ = L(-16/3*(-1+y+2x*(7-6y+x*(-21+16x+12y))), -8/3*x*(-2+4x)*(-1+4x))  
    ∇ϕ₇ = L(16/3*(1+12x*(-1+2x))*y, 8/3*x*(-2+4x)*(-1+4x))  
    ∇ϕ₈ = L(4*(-1+8x)*y*(-1+4y), 4x*(-1+4x)*(-1+8y))  
    ∇ϕ₉ = L(8/3*y*(-2+4y)*(-1+4y), 16/3*x*(1+12y*(-1+2y)))
    ∇ϕ₁₀ = L(-8/3*y*(-2+4y)*(-1+4y), -16/3*(-1+x+12*x*y*(-1+2y)+2*y*(7+y*(-21+16y))))
    ∇ϕ₁₁ = L(4*y*(-1+4y)*(-7+8x+8y), 4*(-1+x+2y)*(3+32*(-1+y)*y+4x*(-1+8y)))
    ∇ϕ₁₂ = L(-16/3*y*(13+24x^2+12y*(-3+2*y)+12x*(-3+4y)), -16/3*(-3+8x^3+6x^2*(-3+8y)+x*(13+72*(-1+y)*y)+2y*(13+y*(-27+16y))))
    ∇ϕ₁₃ = L(32y*(3+2x*(-7+6x)-7y+16*x*y+4y^2), 32x*(3+4x^2+2y*(-7+6y)+x*(-7+16y)))   
    ∇ϕ₁₄ = L(-32y*(1-y+2x*(-5+6x+4y)), -32x*(-1+4x)*(-1+x+2y))  
    ∇ϕ₁₅ = L(-32y*(-1+2x+y)*(-1+4y), -32x*(1+2y*(-5+6y)+x*(-1+8y)))
    v = vals(el)
    return v[1]*∇ϕ₁ + v[2]*∇ϕ₂ + v[3]*∇ϕ₃ +
           v[4]*∇ϕ₄ + v[5]*∇ϕ₅ + v[6]*∇ϕ₆ +
           v[7]*∇ϕ₇ + v[8]*∇ϕ₈ + v[9]*∇ϕ₉ +
           v[10]*∇ϕ₁₀ + v[11]*∇ϕ₁₁ + v[12]*∇ϕ₁₂ + 
           v[13]*∇ϕ₁₃ + v[14]*∇ϕ₁₄ + v[15]*∇ϕ₁₅
end

# P9 for ReferenceTriangle
# source: Silvester PP, Ferrari RL, Finite elements for electrical engineers (1990).
@fastmath function _triangle_rfunc(x, n, m)
    r = 1.0
    for k in 0:m-1
        r *= (x - k/n)/(m/n - k/n) 
    end
    return r
end
@fastmath function _triangle_afunc(x, y, n, i, j, k)
    a = _triangle_rfunc(x, n, i)
    a *= _triangle_rfunc(y, n, j)
    a *= _triangle_rfunc(1-x-y, n, k)
    return a
end

@fastmath function (el::LagrangeTriangle{55})(u)
    @assert u ∈ domain(el)
    x = u[1]
    y = u[2]
    n = 9    # order
    afunc = (i,j,k) -> _triangle_afunc(x,y,n,i,j,k)
    ϕ = (afunc(0,0,9),afunc(9,0,0),afunc(0,9,0),afunc(1,0,8),afunc(2,0,7),
         afunc(3,0,6),afunc(4,0,5),afunc(5,0,4),afunc(6,0,3),afunc(7,0,2),
         afunc(8,0,1),afunc(8,1,0),afunc(7,2,0),afunc(6,3,0),afunc(5,4,0),
         afunc(4,5,0),afunc(3,6,0),afunc(2,7,0),afunc(1,8,0),afunc(0,8,1),
         afunc(0,7,2),afunc(0,6,3),afunc(0,5,4),afunc(0,4,5),afunc(0,3,6),
         afunc(0,2,7),afunc(0,1,8),afunc(1,1,7),afunc(7,1,1),afunc(1,7,1),
         afunc(2,1,6),afunc(3,1,5),afunc(4,1,4),afunc(5,1,3),afunc(6,1,2),
         afunc(6,2,1),afunc(5,3,1),afunc(4,4,1),afunc(3,5,1),afunc(2,6,1),
         afunc(1,6,2),afunc(1,5,3),afunc(1,4,4),afunc(1,3,5),afunc(1,2,6),
         afunc(2,2,5),afunc(5,2,2),afunc(2,5,2),afunc(3,2,4),afunc(4,2,3),
         afunc(4,3,2),afunc(3,4,2),afunc(2,4,3),afunc(2,3,4),afunc(3,3,3))
    v = vals(el)
    @assert 55 == length(v) == length(ϕ)
    result = first(v) * first(ϕ)
    for i in 2:55
        result += v[i]*ϕ[i]
    end
    return result
end

function jacobian(el::LagrangeTriangle{55}, u)
    @assert u ∈ domain(el)
    return ForwardDiff.jacobian(el,u)
end
