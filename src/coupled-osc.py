from scipy.integrate import odeint
import numpy as np

A = 1
B = 1


def brusselator(c, t):
    x = c[0]
    y = c[1]

    dxdt = A - (B-1)*x
    dydt = B*x - x*x*y

    return [dxdt, dydt]


def brusselator_gradient(c, t):
    x = c[0]
    y = c[1]

    dxdx = -(B-1)
    dxdy = 0
    dydx = B - 2*x*y
    dydy = -x*x

    return [[dxdx, dxdy], [dydx, dydy]]

t = np.arange(0, 4.0, 0.01)

c0 = [1.0, 1.0]
c = odeint(brusselator, c0, t, Dfun=brusselator_gradient)
print(c)
