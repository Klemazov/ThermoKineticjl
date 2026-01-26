#Kintic ODE Function 
# it should be a storage for creating K(T(t))*f(α)

struct KineticODEFuncion
    K
    A
    Ea
    f
    T
    params
end
"""
    (k::KineticODEFuncion)(t,α)

TBW
"""
function (k::KineticODEFuncion)(t,α)
    #try to create dαdt = K(T(t))*f(α)
    K = k.T
    A = k.A
    Ea = k.Ea
    f = k.f
    T = k.T
    params = k.params
    return K(A, Ea,T(t))*f(α, params...)
end


