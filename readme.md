# CasADi symbolic splines
The library helps to create symbolic expressions from B-Splines.

# Example:
```python
from casadi_bspline.bspline import get_spline_symfun
from scipy.interpolate import make_interp_spline
import numpy as np
import casadi as ca

x = np.linspace(-1, 2, 100)
y = np.sin(x)
sp = make_interp_spline(x, y, k=5)
spline_fun = get_spline_symfun(sp, ca.SX)
spline_fun_deriv = spline_fun.jacobian()

x0 = 0.729364
print(f'spline value at {x0} is {spline_fun(x0)}')
print(f'spline derivative value at {x0} is {spline_fun_deriv(x0, 1)}')
```

# Test
python -m unittest -v src/test.py
