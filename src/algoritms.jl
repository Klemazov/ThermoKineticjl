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

res_rk = rk4(f, (0.,1.0), 1.0, 10)
res_euler = euler(f, (0.0,1.0), 1.0, 10)
res_midpoint = midpoint(f, (0.,1.0), 1.0, 10)

function plot_res()
    plot(range(0,1,100), exp.(range(0,1, 100)))
    scatter!(res_rk, label = "rk4")
    scatter!(res_euler, label = "euler")
    scatter!(res_midpoint, label = "mp")
end 


#= passing function as arguments

1) direct function call
    function g_inline(f) 
        a = rand()
        return f(a)
    end
2) indirect function call
 function g_indirect(f)
        a = rand()
        apply_func(f,a)
  end
  apply_func(f,a) = f(a)

3) closure call
  function g_cloure(f)
    a = rand()
    return apply_func_closure(a)
  end
  apply_func_closure(a) = f(a)


  1) - fast 
  2) - slow because internal function cnnot be specialized on f type at copilation time
  3) - fast because apply_func_closure is inlined 
=#
# begin
# using BenchmarkTools
# function make_grid(u₀, tspan, n)
#     t0, t_end = tspan
#     h = (t_end-t0)/n
#     t = [t0+i*h for i in 0:n]
#     u = fill(u₀..., n+1)
#     return (t, h, u)
# end

# # indirect function call
# function rk4_indirect(f, tspan, u₀, n)
#     (t, h, u) = make_grid(u₀, tspan, n)
#     for i in 1:n
#        @inbounds u[i+1] = u[i] +  rk4_iterate_indirect(f,t[i],h,u[i]) # passing function as internal function argument
#     end
#     return t,u
# end
# function rk4_iterate_indirect(f,t,h,u) # this function doesnt know the 
#     # type of f at compile time! 
#         k₁ = f(t, u)
#         k₂ = f(t+0.5h, u+0.5h*k₁)
#         k₃ = f(t, u+0.5h*k₂)
#         k₄ = f(t, u+h*k₃)
#         return (h/6)*(k₁+2k₂+2k₃+k₄)
# end

# # closure function call
# function rk4_closure(f, tspan, u₀, n)
#     (t, h, u) = make_grid(u₀, tspan, n)
#     for i in 1:n
#         @inbounds u[i+1] = u[i] +  rk4_iterate_closure(t[i],h,u[i])
#     end
#     return t,u
# end
# function rk4_iterate_closure(t,h,u)
#     # uses external function f specification
#         k₁ = f(t, u)
#         k₂ = f(t+0.5h, u+0.5h*k₁)
#         k₃ = f(t, u+0.5h*k₂)
#         k₄ = f(t, u+h*k₃)
#         return (h/6)*(k₁+2k₂+2k₃+k₄)
# end


# # closure type wrapper
# function rk4_closure_wrapper_type(f, tspan, u₀, n)
#     (t, h, u) = make_grid(u₀, tspan, n)
#     c_w = call_wrapper_wrapper() # closure wrapper wrapps around runtime type of f
#     for i in 1:n
#         @inbounds u[i+1] = u[i] +  rk4_iterate_closure_wrapper_type(c_w,t[i],h,u[i]) # passing function as internal function argument
#     end
#     return t,u
# end
# struct call_wrapper{F}
#     func::F
# end
# call_wrapper_wrapper() = call_wrapper(f) # this function must be called within the method 
# # to infer the current function type 
# function rk4_iterate_closure_wrapper_type(c_w::call_wrapper{F},t,h,u) where F
#         k₁ = c_w.func(t, u)
#         k₂ = c_w.func(t+0.5h, u+0.5h*k₁)
#         k₃ = c_w.func(t, u+0.5h*k₂)
#         k₄ = c_w.func(t, u+h*k₃)
#         return (h/6)*(k₁+2k₂+2k₃+k₄)
# end


# print("Direct call all inline:")
# @btime rk4(f, (0.,1.0), 1.0, 10)

# print("Indirect call with passing function directly as an argument: ")
# @btime rk4_indirect(f, (0.,1.0), 1.0, 10)

# print("Indirect call using closure pattern: ")
# @btime rk4_closure(f, (0.,1.0), 1.0, 10)

# print("Indirect call with closure wrapping using struct ")
# @btime rk4_closure_wrapper_type(f, (0.,1.0), 1.0, 10)

# end;