
abstract type ReactionRateModel end


function K(reaction_rate::ReactionRateModel)
    A = reaction_rate.A
    Ea = reaction_rate.Ea
    T = reaction_rate.T
    return A*exp(-Ea/(R*T))
end

struct ReactionRate{T} <: ReactionRateModel
    A::T
    Ea::T
    T::T
end

Ka = ReactionRate(1e5,1e5,1000.0)

K(Ka)
# K(A, Ea,T) = A*exp(-Ea/(R*T))
# K(A, Ea,T, ::Type{ReactionRateA}) = K(A, Ea,T)

# K(B,Ea,T,::Type{ReactionRateB}) = B*√T*exp(-Ea/(R*T))
# K(C, Ea, T, ::Type{ReactionRateC}) = C*T*exp(-Ea/(R*T))

