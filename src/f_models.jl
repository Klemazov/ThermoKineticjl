
f1(α) = 1-α

f2(α) = (1-α)^2

fn(α,n) = (1-α)^n

pt(α, n, m) = α^n*(1-α)^m

ks(α,n,m,Kcat) = fn(α, n)*(1+Kcat*α^m)



#define base structures

abstract type FunctionModel end

mutable struct Fn<:FunctionModel
    params::NTuple{1,Float64}
    name::NTuple{1,Symbol}

    Fn(n=1.0) = begin
         new(Tuple(n,), (:n,))
    end
end

mutable struct F1<:FunctionModel
    name::NTuple{1,Symbol}
    F1() = begin
         new((:one,))
    end
end

mutable struct F2<:FunctionModel
    name::NTuple{1,Symbol}
    F2() = begin
         new((:two,))
    end
end

mutable struct Pt<:FunctionModel
    params::NTuple{2,Float64}
    name::NTuple{2,Symbol}

    Pt(n=1.0, m=1.0) = begin
         new((n,m), (:n,:m))
    end
end

mutable struct Ks<:FunctionModel
    params::NTuple{3,Float64}
    name::NTuple{3,Symbol}

    Ks(n=1.0, m=1.0, Kcat = 1.0) = begin
         new((n,m, Kcat), (:n,:m, :Kcat))
    end
end

#create functors
 for (model_name,f_name) in zip((:F1,:F2,),
                            (:f1, :f2,))
        @eval (p::$model_name)(α::T) where T= $f_name(α)::T
 end



 for (model_name,f_name) in zip((:Fn,:Pt,:Ks),
                            (:fn, :pt, :ks))
        @eval (p::$model_name)(α::T) where T= $f_name(α,getfield(p,:params)...)::T
 end

