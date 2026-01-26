module ThermalKineticjl
include("algoritms.jl")
include("f_models.jl")
include("reaction_constant_model.jl")
include("temperature_functions.jl")
include("ThermalKineticTypes.jl")
# Write your package code here.


const R = 8.314

exemplar = KineticODEFuncion(K,1,1,f1,T,[1])


display(exemplar(1,0.1))
end
