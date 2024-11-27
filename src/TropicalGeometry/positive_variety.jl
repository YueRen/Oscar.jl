@doc raw"""
    positive_tropical_variety(I::MPolyIdeal,nu::TropicalSemiringMap)

Return the positive tropical variety of `I` as a `PolyhedralComplex` as per the definition in [SW05](@cite).  Assumes that `I` is generated by affine linear polynomials either
(a) defined over the rational numbers and that `nu` encodes the trivial valuation,
(b) defined over the rational function field over the rational numbers and that `nu` encodes the t-adic valuation.

# Examples
```jldoctest
julia> K,t = rational_function_field(QQ,"t")
(Rational function field over QQ, t)

julia> C = matrix(K,[[-3*t,1*t,-1*t,-2*t,2*t],[-1*t,1*t,-1*t,-1*t,1*t]])
[-3*t   t   -t   -2*t   2*t]
[  -t   t   -t     -t     t]

julia> R,x = polynomial_ring(K,ncols(C))
(Multivariate polynomial ring in 5 variables over K, AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.RationalFunctionFieldElem{QQFieldElem, QQPolyRingElem}}[x1, x2, x3, x4, x5])

julia> nu = tropical_semiring_map(K,t)
Map into Min tropical semiring encoding the t-adic valuation on Rational function field over QQ

julia> I = ideal(C*gens(R))
Ideal generated by
  -3*t*x1 + t*x2 - t*x3 - 2*t*x4 + 2*t*x5
  -t*x1 + t*x2 - t*x3 - t*x4 + t*x5

julia> TropPlusI = positive_tropical_variety(I,nu)
Min tropical variety

```
"""
function positive_tropical_variety(I::MPolyIdeal,nu::TropicalSemiringMap)
    @req all(isequal(1),total_degree.(gens(I))) "generators of input ideal not affine linear"

    # Construct the tropicalization of I
    TropL = tropical_linear_space(I,nu)

    # find maximal polyhedra belonging to the positive part
    positivePolyhedra = Polyhedron{QQFieldElem}[]
    for sigma in maximal_polyhedra(TropL)
        # compute a reduce Groebner basis with respect to a relative interior point
        w = relative_interior_point(sigma)
        inI = initial(I,nu,w)
        G = groebner_basis(inI; complete_reduction=true)

        # the Groebner basis is binomial, check binomials have alternating signs
        @req all(isequal(2), length.(G)) "initial ideal is not binomial"
        if all(isequal(-1),[prod([sign(c) for c in coefficients(g)]) for g in G])
            push!(positivePolyhedra,sigma)
        end
    end

    if isempty(positivePolyhedra)
        # if there are no positive polyhedra,
        # return empty polyhedral complex in the correct ambient dimension
        return polyhedral_complex(IncidenceMatrix(zeros(Int,0,0)),zero_matrix(QQ,0,ambient_dim(TropL)))
    end

    Sigma = polyhedral_complex(positivePolyhedra)
    mult = ones(ZZRingElem, n_maximal_polyhedra(Sigma))
    minOrMax = convention(nu)
    return tropical_variety(Sigma,mult,minOrMax)
end
