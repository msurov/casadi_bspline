import unittest
from bspline import (
    get_spline_symexpr,
    get_spline_sx,
    get_spline_mx,
    get_spline_symfun
)
from scipy.interpolate import make_interp_spline
import numpy as np
import casadi as ca



class TestClass(unittest.TestCase):
        
    def __spline_expr(self, symtype):
        x = np.linspace(-1, 2, 100)
        y = np.sin(x)
        sp = make_interp_spline(x, y, k=5)
        xsym = symtype.sym('x')
        splineexpr = get_spline_symexpr(sp, xsym)
        splineexpr_deriv = ca.jacobian(splineexpr, xsym)

        xmin = np.min(x)
        xmax = np.max(x)
        x_test = np.random.rand(100) * (xmax - xmin) + xmin

        for x0 in x_test:
            val = ca.substitute(splineexpr, xsym, x0)
            val = float(ca.evalf(val))
            self.assertAlmostEqual(val, sp(x0))

            val = ca.substitute(splineexpr_deriv, xsym, x0)
            val = float(ca.evalf(val))
            self.assertAlmostEqual(val, sp(x0, 1))

    def __spline_fun(self, symtype):
        x = np.linspace(-1, 2, 100)
        y = np.sin(x)
        sp = make_interp_spline(x, y, k=5)
        spline_fun = get_spline_symfun(sp, symtype)
        spline_fun_deriv = spline_fun.jacobian()

        xmin = np.min(x)
        xmax = np.max(x)
        x_test = np.random.rand(100) * (xmax - xmin) + xmin

        for x0 in x_test:
            val = spline_fun(x0)
            val = float(ca.evalf(val))
            self.assertAlmostEqual(val, sp(x0))

            val = spline_fun_deriv(x0, 1)
            val = float(ca.evalf(val))
            self.assertAlmostEqual(val, sp(x0, 1))

    def __extrapolate(self, symtype):
        x = np.linspace(-1, 2, 100)
        y = np.sin(x)
        sp = make_interp_spline(x, y, k=5)
        spline_fun = get_spline_symfun(sp, symtype)
        spline_fun_deriv = spline_fun.jacobian()

        xmin = np.min(x)
        xmax = np.max(x)

        x_test = [
            xmin + (xmax - xmin) * 1.521,
            xmin + (xmax - xmin) * -0.22456
        ]

        for x0 in x_test:
            val = spline_fun(x0)
            val = float(ca.evalf(val))
            self.assertAlmostEqual(val, sp(x0))

            val = spline_fun_deriv(x0, 1)
            val = float(ca.evalf(val))
            self.assertAlmostEqual(val, sp(x0, 1))
            
    def test_extrapolate_sx(self):
        self.__extrapolate(ca.SX)

    def test_extrapolate_mx(self):
        self.__extrapolate(ca.MX)

    def test_spline_expr_sx(self):
        self.__spline_expr(ca.SX)

    def test_spline_expr_mx(self):
        self.__spline_expr(ca.MX)

    def test_spline_fun_sx(self):
        self.__spline_fun(ca.SX)

    def test_spline_fun_mx(self):
        self.__spline_fun(ca.MX)

if __name__ == '__main__':
    unittest.main()
