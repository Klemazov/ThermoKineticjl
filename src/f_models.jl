
"""
    first order reaction
    f1(α)

TBW
"""
f1(α) = 1-α

"""
    second order reaction
    f2(α)

TBW
"""
f2(α) = (1-α)^2

"""
    n-th order reaction
    fn(α,n)

TBW
"""
fn(α,n) = (1-α)^n


"""
    Prout-Tompkins reaction
    pt(α)

TBW
"""
pt(α) = α(1-α)


"""
    expanded Prout-Tompkins reaction
    pt_expand(α, n, m)

TBW
"""
pt_expand(α, n, m) = α^n*(1-α)^m


"""
    Kama-Sourour auto-catalysis
    ks(α,n,m,Kcat)

TBW
"""
ks(α,n,m,Kcat) = fn(α, n)*(1+Kcat*α^m)


