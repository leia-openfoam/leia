import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def blasius_ode(eta, y):
    # y = [f, f', f'']
    return [y[1], y[2], -0.5*y[0]*y[2]]

def shoot(s, eta_max=10.0):
    # integrate with initial slope y2(0)=s, return f'(eta_max)-1
    y0 = [0.0, 0.0, s]
    sol = solve_ivp(blasius_ode, [0, eta_max], y0, dense_output=True, max_step=0.1)
    return sol.y[1, -1] - 1.0

# find the correct second derivative at the wall
s_star = brentq(shoot, 0.3, 0.4)     # bracket around ~0.332
print(f"Found f''(0) ≈ {s_star:.6f}")

# integrate with the correct shooting parameter
eta = np.linspace(0, 10, 200)
y0 = [0.0, 0.0, s_star]
sol = solve_ivp(blasius_ode, [0, 10], y0, t_eval=eta, max_step=0.05)

plt.figure(figsize=(6,4))
plt.plot(sol.t, sol.y[1], lw=2)
plt.xlabel(r'$\eta$')
plt.ylabel(r'$f\'(\eta)\ \equiv\ u/U_\infty$')
plt.title('Blasius Velocity Profile')
plt.grid(True)
plt.ylim(0, 1.05)
plt.show()