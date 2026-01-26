
struct ReactionRateA end

struct ReactionRateB end

struct ReactionRateC end

K(A, Ea,T) = A*exp(-Ea/(R*T))
K(A, Ea,T, ::Type{ReactionRateA}) = K(A, Ea,T)

K(B,Ea,T,::Type{ReactionRateB}) = B*√T*exp(-Ea/(R*T))
K(C, Ea, T, ::Type{ReactionRateC}) = C*T*exp(-Ea/(R*T))

