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

f(t, u) = u;

res_rk = rk4(f, (0.,10.0), 1.0, 100);
res_euler = euler(f, (0.0,10.0), 1.0, 100);
res_midpoint = midpoint(f, (0.,10.0), 1.0, 100);

function plot_res()
    plot(range(1,10,1000), exp.(range(1,10,1000)))
    scatter!(res_rk, label = "rk4")
    scatter!(res_euler, label = "mp")
    scatter!(res_midpoint, label = "euler")
end 