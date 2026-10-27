#!/usr/bin/env python3
"""Two measures of the normal error, symbolically (Appendix A of the manuscript).

A recovered unit normal n~ that deviates from the exact unit normal n by the angle theta
has |n - n~| = 2 sin(theta/2) and 1 - n . n~ = 1 - cos(theta). The series in theta show
that the first measure is linear in the angle and the second quadratic, so a normal whose
angle error is O(h) (first order) appears as O(h^2) in the second measure, in L1 and L2
norms alike (the norms of theta^2 over the cells are second order when theta is first).
Run: python3 normal_error_measures_sympy.py
"""
import sympy as sp

theta, h, c = sp.symbols("theta h c", positive=True)
n = sp.Matrix([0, 0, 1])
n_tilde = sp.Matrix([sp.sin(theta), 0, sp.cos(theta)])  # the exact normal turned by theta

d_vec = sp.simplify(sp.sqrt((n - n_tilde).dot(n - n_tilde)))
d_dot = sp.simplify(1 - n.dot(n_tilde))
print("|n - n~|      =", sp.simplify(d_vec), " = ", sp.series(d_vec, theta, 0, 6).removeO())
print("1 - n . n~    =", d_dot, " = ", sp.series(d_dot, theta, 0, 6).removeO())

# a first-order normal: theta = c h
order_vec = sp.limit(sp.log(d_vec.subs(theta, c * h)) / sp.log(h), h, 0)
order_dot = sp.limit(sp.log(d_dot.subs(theta, c * h)) / sp.log(h), h, 0)
print("order in h of |n - n~| for theta = c h:", order_vec)
print("order in h of 1 - n.n~ for theta = c h:", order_dot)
# the ratio of the two measures at a given angle
print("(1 - n.n~) / |n - n~|^2 =", sp.simplify(d_dot / d_vec**2))
