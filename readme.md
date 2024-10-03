# CasADi symbolic splines
The library helps to create symbolic expressions from B-Splines.

# Install
Install from github
```bash
python3 -m pip install git+https://github.com/msurov/casadi_bspline.git
```
Clone and run the tests
```bash
git clone https://github.com/msurov/casadi_bspline.git
cd casadi_bspline
python3 -m venv .venv
python3 -m pip install .
python3 -m unittest src/test.py
```

# Example
```python
from bspline import get_spline_symfun
from scipy.interpolate import make_interp_spline
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 2*np.pi, 100)
y = np.sin(x)
dy = np.cos(x)
sp = make_interp_spline(x, y, k=5, bc_type='periodic')
spline_fun = get_spline_symfun(sp, ca.SX)
spline_fun_deriv = spline_fun.jacobian()
deriv_vals = spline_fun_deriv(x, 1)
plt.plot(x, deriv_vals, lw=2, color='b')
plt.plot(x, dy, '--', lw=2, color='red', alpha=1)
plt.grid(True)
```
