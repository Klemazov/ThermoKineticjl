begin
struct ParamsFunctionHandler{F,P}
    f::F
    params::P
end
"""

    (handler::ParamsFunctionHandler)(u,t)

TBW
"""
function (handler::ParamsFunctionHandler)(t,u)
    handler.f(t,u,handler.params...)
end


struct IVPFunction
    f
    params
end

struct IVP{T,F}
    tspan::NTuple{2,T}
    u₀::T
    f::F
end

"""
    euler(f, tspan, u₀, n)


TBW
"""
function euler(f, tspan, u₀, n)
    t0, t_end = tspan;
    h = (t_end-t0)/n;
    t = [t0+i*h for i in 0:n]
    u = fill(u₀..., n+1)

    for i in 1:n
        u[i+1] = u[i]+h*f(t[i],u[i])
    end
    return t,u
end

"""
    midpoint(f, tspan, u₀, n)


TBW
"""
function midpoint(f, tspan, u₀, n)
    t0, t_end = tspan;
    h = (t_end-t0)/n;
    t = [t0+i*h for i in 0:n]
    u = fill(u₀..., n+1)

    for i in 1:n
        k₁ = f(t[i], u[i])
        uₘ = u[i] + 0.5*h*k₁;
        tₘ = t[i]+0.5h;
        k₂ = f(tₘ, uₘ);
        u[i+1] = u[i] + h* k₂
    end
    return t,u
end

"""
    rk4(f, tspan, u₀, n)


TBW
"""
function rk4(f, tspan, u₀, n)
    t0, t_end = tspan;
    h = (t_end-t0)/n;
    t = [t0+i*h for i in 0:n]
    u = fill(u₀..., n+1)

    for i in 1:n
        k₁ = f(t[i], u[i])
        k₂ = f(t[i]+0.5h, u[i]+0.5h*k₁)
        k₃ = f(t[i], u[i]+0.5h*k₂)
        k₄ = f(t[i], u[i]+h*k₃)
        u[i+1] = u[i] + (h/6)*(k₁+2k₂+2k₃+k₄)
    end
    return t,u
end

function f_t(t,u)
    T = 9t
    return exp(-1/T)*(1-u)^5
end

f(t, u) = 0.2exp(-20/t)*(1-u)^8+0.8exp(-90/t)*(1-u)^2;
f1(t,u) = exp(-1/t)*(1-u)
fcat(t,u) = exp(-1/t)*0.2(1-u)*u

res_f = rk4(f, (0.1,100.0), 0.001, 1000);
res_f1=rk4(f1, (0.1,100.0), 0.001, 1000);
res_cat = res_rk_f1=rk4(fcat, (0.1,100.0), 0.001, 1000);
res_f_t = rk4(f_t, (0.1,100.0), 0.001, 1000);

#plot(res_f)
#plot!(res_f1)
#plot!(res_cat)
plot(res_f_t)
end