from casadi import SX, MX, if_else, DM, Function, substitute, reshape
from scipy.interpolate import BSpline
import numpy as np


def __select_value(x, t, c):
    """
        Binary search. Returns symbolic expression which selects c[j] if t[j] <= x < t[j+1]

        Keyword arguments:
            x -- symbolic argument to search for
            t -- knots
            c -- control points
    """
    N = len(t)
    if N == 1:
        return c[0]
    return if_else(
        x >= t[N//2],
        __select_value(x, t[N//2:], c[N//2:]),
        __select_value(x, t[:N//2], c[:N//2]),
    )

def __bisect_expr(arg, knots : DM):
    '''
        Get symbolic expression for the bisect algorithm

        Keyword arguments:
            arg -- symbolic argument to search for
            knots -- array of spline knots
    '''
    N = knots.shape[0]
    if N == 1:
        return 0
    m = N // 2
    return if_else(
        arg >= knots[m],
        m + __bisect_expr(arg, knots[m:]),
        __bisect_expr(arg, knots[:m])
    )

def __compute_basis(x, i, k, t):
    '''
        The method computes of basis functions 
            B_{i-k,k}, B_{i-k+1,k}, ..., B_{i,k}

        Keyword arguments:
            x -- symbolic variable, it is assumed it belongs to the half-interval [t_{i} .. t_{i+1}) 
            i -- the index of the knots interval  [t_{i} .. t_{i+1})  containing x
            k -- spline degree
            t -- spline knots
    '''
    sym_type = type(x)
    D = sym_type.zeros(k + 1, k + 1)
    D[0,k] = 1

    for r in range(1, k + 1):
        D[r, k-r] = (t[i+1] - x) / (t[i+1] - t[i-r+1]) * D[r-1, k-r+1]
        for a in range(1, r):
            D[r, k-a] = \
                (x - t[i-a]) / (t[i-a+r] - t[i-a]) * D[r-1, k-a] + \
                (t[i-a+r+1] - x) / (t[i-a+r+1] - t[i-a+1]) * D[r-1, k-a+1]
        D[r, k] = (x - t[i]) / (t[i+r] - t[i]) * D[r-1, k]

    Bk = D[k,:].T
    return Bk


def get_spline_sx(sp : BSpline, x : SX):
    '''
        Makes symbolic expression for the given BSpline
        TODO: add extrapolation type arguments
    '''
    t = sp.t
    assert isinstance(x, SX), 'Support only SX and MX types for the spline argument'
    sym_type = SX

    N,*dshape = sp.c.shape
    if dshape == []:
        dshape = [1,1]
    elif len(dshape) == 1:
        dshape = dshape + [1,]
    total = np.prod(dshape, dtype=int)
    c = np.reshape(sp.c, (N, total))
    k = sp.k
    result = 0

    for i in range(k, N):
        if i == k:
            nonzero = x < t[i + 1]
        elif i == N - 1:
            nonzero = x >= t[i]
        else:
            nonzero = (x >= t[i]) * (x < t[i + 1])
        result += nonzero * __compute_basis(x, i, k, t).T @ c[i-k:i+1,:]

    return result


def get_spline_mx(sp : BSpline, x : MX):
    '''
        Makes symbolic expression for the given BSpline
        TODO: add extrapolation type arguments
    '''
    t = sp.t
    assert isinstance(x, MX), 'Support only SX and MX types for the spline argument'
    sym_type = MX

    N,*dshape = sp.c.shape
    if dshape == []:
        dshape = [1,1]
    elif len(dshape) == 1:
        dshape = dshape + [1,]
    total = np.prod(dshape, dtype=int)
    c = np.reshape(sp.c, (N, total))
    k = sp.k
    values = sym_type.zeros(N-k, total)

    for i in range(k, N):
        values[i - k,:] = __compute_basis(x, i, k, t).T @ c[i-k:i+1,:]

    idx = __bisect_expr(x, t[k:N])
    expr = values[idx,:]
    expr = reshape(expr, *dshape[::-1]).T
    return expr

def get_spline_symexpr(sp : BSpline, x : SX | MX):
    if isinstance(x, SX):
        return get_spline_sx(sp, x)
    elif isinstance(x, MX):
        return get_spline_mx(sp, x)
    else:
        assert False, 'x expected to be of SX|MX type'

def get_spline_symfun(sp : BSpline, symtype = SX):
    if symtype == SX:
        x = SX.sym('dummy')
        expr = get_spline_sx(sp, x)
    elif symtype == MX:
        x = MX.sym('dummy')
        expr = get_spline_mx(sp, x)
    else:
        assert False, 'x expected to be of SX|MX type'

    return Function('Spline', [x], [expr])
