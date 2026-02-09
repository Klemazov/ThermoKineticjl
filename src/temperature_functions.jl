
abstract type TemperatureModel end

struct TConstant<:TemperatureModel end
struct TConstantRate<:TemperatureModel end
struct TPolynomial<: TemperatureModel end




T(t) = t

T(t,β) = β*T

#TODO
#T(t) = spline(t)  