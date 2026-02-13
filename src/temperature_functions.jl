
abstract type TemperatureModel end

struct Temperature{F}<:TemperatureModel
    f::F
end

(p::Temperature)() = p.f

