
abstract type KModel end

struct  KineticRate{P}<:KModel 
    Ea::P
    A::P
end

function (k::KineticRate)(T)
     A = k.A
     Ea = k.Ea
    return A*exp(-Ea/(R*T))
end

