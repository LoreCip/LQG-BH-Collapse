import time
import numpy as np

import matplotlib.pyplot as plt

from scipy.integrate import cumulative_trapezoid
from numba import njit, vectorize

# @njit
def f(u, x):
    return 0.5 * x*x*x * np.sin(u / x/x)**2 

# @njit
def fp(u, x):
    return 0.5 * x * np.sin(2 * u /x/x)

def initial_data(x, m, r0, h):

    # rho = 3 * m * (1 - heaviside(xs-r0)) / (4 * np.pi * r0**3)
    # rho = 3 * m / (4*np.pi * r0**3) * (np.pi/2 - np.arctan(x - r0))/2
    
    # M = 4*np.pi * cumulative_trapezoid(rho * x*x, x, initial=0)
    # print(f'Mass computed at infinity = {np.round(M[-1],3)}')

    M = m * x**3 / r0**3 * heaviside(r0 - x) + m * heaviside(x - r0)

    B0 = - 0.5 * x**2 * np.arccos(1 - 4 * M / x**3)
    B0[0] = 0
    
    return B0

# @vectorize
# def heaviside(x):
#     if x < 0:
#         return 0
#     elif x > 0:
#         return 1
#     else:
#         return 0.5

def heaviside(x):
    out = np.zeros_like(x)
    for i, xx in enumerate(x):
        if xx > 0:
            out[i] = 1
        elif xx < 0:
            out[i] = 0
        else:
            out[i] = 0.5
    return out

# @njit
def TDV_RK(r, u_prev, xs, h, dt, nghost):
    if r == 2:
        return TDV_RK3(u_prev, xs, h, dt, nghost)
    else:
        raise Exception("Order not implemented! r = 2")

# @njit
@profile
def L(r, j, u, fup, xs, h, nghost):
    """
    Equation 2.7b

    u  --> array of sol
    xs --> array of grid
    
    x_jm = x_{j - 1/2}
    x_jM = x_{j + 1/2}

    R(r, j  , x_jM, u, fup, xs, h) = f(x_{j + 1/2})
    """
    x_jM = xs[j] + h / 2
    x_jm = xs[j] - h / 2
    if j == nghost:
        o1 = flux(R1, R2, x_jM) 
        o2 = 0
    else:
        o1 = flux(R(r, j  , x_jM, u, fup, xs, h), R(r, j+1, x_jM, u, fup, xs, h), x_jM) 
        o2 = flux(R(r, j-1, x_jm, u, fup, xs, h), R(r, j, x_jm, u, fup, xs, h), x_jm)
    return - ( o1 - o2 ) / h

# @njit
def R(r, j, x_pt, u, fup, xs, h):
    if r == 2:
        return R2(j, x_pt, u, fup, xs, h)
   
# @njit
def alpha(r, j, i, u, fup, x, h):
    if r == 2:
        return alpha_r2(j, i, u, fup, x, h)

# @njit
def interpolants(r, j, x, u, xs, h):
    if r == 2:
        return interp_r2(j, x, u, xs, h)


# @njit
def SI(r, j, u):
    if r == 2:
        return SI_r2(j, u)

# @njit
# def flux(a, b, x):
#     """
#     Roe flux function
#     Equation 2.6
#     """
#     u_dec = np.linspace(min(a,b), max(a,b))
#     der = fp(u_dec, x)
#     if np.all(der >= 0):
#         return f(a, x)
#     elif np.all(der < 0):
#         return f(b, x)
#     else:
#         beta = np.max(np.abs(der)) 
#         return 0.5* (f(a, x) + f(b, x) - beta * (b - a))

# @njit
# def flux(a, b, x):
#     """
#     Godunov flux function
#     Equation 2.5
#     """
#     u_dec = np.linspace(min(a,b), max(a,b), 100)
#     ff = f(u_dec, x)
#     if a <= b:
#         return np.min(ff)
#     elif a > b:
#         return np.max(ff)

# @njit
@profile
def flux(a, b, x):
    ul = a/x/x
    ur = b/x/x
    FL = 0.5 * x**3 * np.sin(ul)**2 
    FR = 0.5 * x**3 * np.sin(ur)**2 
    if ul <= ur:
        return min(FL, FR)
    elif ul > ur:
        if ((ur > -np.pi/2) or (ul < -np.pi/2)):
            return max(FL, FR)
        else:
            return 0.5 * x**3

# @njit
def BC(arr, nghost):
    for i in range(nghost):
        arr[i] = arr[nghost + 1 - i]
    # arr[:nghost] = arr[nghost:2*nghost]
    arr[-nghost:] = arr[-nghost - 1]
    return arr

@profile
# @njit
def TDV_RK3(u_prev, xs, h, dt, nghost):
    """
    ## TO OPTIMIZE

    1**) Euler step t      --> t +   dt
    2**) Euler step t + dt --> t + 2*dt
    
    3**) Weighted mean btw t and t+2*dt to find t + dt/2
    
    4**) Euler step t + dt/2 -- > t + 3*dt/2
    
    5**) Weighted mean btw t and t+3*dt/2 to find  t + dt
    """
    rr = np.zeros(len(u_prev) - 2*nghost)
    u_n1    = np.zeros_like(u_prev)
    u_n2    = np.zeros_like(u_prev)
    u_n12   = np.zeros_like(u_prev)
    u_n32   = np.zeros_like(u_prev)
    u_final = np.zeros_like(u_prev)

    fup = fp(u_prev, xs)
    # 1**)
    for j in range(len(rr)-2):
        rr[j] = L(2, j+nghost, u_prev, fup, xs, h, nghost)
    u_n1[nghost] = 0
    u_n1[nghost+1:-nghost] = u_prev[nghost+1:-nghost] + dt * rr[1:]
    u_n1 = BC(u_n1, nghost)
    fup = fp(u_n1, xs)

    # 2**)
    for j in range(len(rr)-2):
        rr[j] = L(2, j+nghost, u_n1, fup, xs, h, nghost)
    u_n2[nghost] = 0
    u_n2[nghost+1:-nghost] = u_n1[nghost+1:-nghost] + dt * rr[1:]
    u_n2 = BC(u_n2, nghost)

    # 3**)
    u_n12[nghost] = 0
    u_n12[nghost+1:-nghost] = 3 * u_prev[nghost+1:-nghost] / 4 + u_n2[nghost+1:-nghost] / 4
    u_n12 = BC(u_n12, nghost)
    fup = fp(u_n12, xs)

    # 4**)
    for j in range(len(rr)-2):
        rr[j] = L(2, j+nghost, u_n12, fup, xs, h, nghost)
    u_n32[nghost] = 0
    u_n32[nghost+1:-nghost] = u_n12[nghost+1:-nghost] + dt * rr[1:]
    u_n32 = BC(u_n32, nghost)

    # 5**)
    u_final[nghost] = 0
    u_final[nghost+1:-nghost] = u_prev[nghost+1:-nghost] / 3 + 2 * u_n32[nghost+1:-nghost] / 3
    u_final = BC(u_final, nghost)

    return u_final

# @njit
@profile
def R2(j, x_pt, u, fup, xs, h):
    """
    Equation 3.16
    """
    aj0 = alpha(2, j, 0, u, fup, x_pt, h)
    aj1 = alpha(2, j, 1, u, fup, x_pt, h)
    w0 = aj0 / (aj0 + aj1)
    w1 = aj1 / (aj0 + aj1)
    out = w0 * interpolants(2, j, x_pt, u, xs, h) + w1 * interpolants(2, j+1, x_pt, u, xs, h)
    return out


#####################
###### OPTIONS FOR 
#####################
# SYMMETRY WRT x_{j}
# For wave left -> right         For wave right -> left
# IS_0 = u[j] - u[j-1]           IS_0 = u[j+1] - u[j]
# IS_1 = u[j+1] - u[j]           IS_1 = u[j] - u[j-1]

# SYMMETRY WRT x_{j + 1/2}
# For wave left -> right         For wave right -> left
# IS_0 = u[j] - u[j-1]           IS_0 = u[j+1] - u[j]
# IS_1 = u[j+1] - u[j]           IS_1 = u[j+2] - u[j+1]

# @njit
def alpha_r2(j, i, u, fup, x, h, eps = 1e-6):
    """
    Equations 3.17a & 3.17b 
    """
    if fup[j] > 0:
        if i == 0:
            # return 1 / 2 / (eps + SI(2, j, u))**2
            return 1 / 2 / (eps + (u[j] - u[j-1])**2)**2
        elif i == 1:
            # return 1 / (eps + SI(2, j+1, u))**2
            return 1 / 1 / (eps + (u[j+1] - u[j])**2)**2
    elif fup[j] <= 0:
        if i == 0:
            # return 1 / (eps + SI(2, j, u))**2
            return 1 / 1 / (eps + (u[j] - u[j+1])**2)**2
        elif i == 1:
            # return 1 / 2 / (eps + SI(2, j+1, u))**2 
            return 1 / 2 / (eps + (u[j+1] - u[j+2])**2)**2

# @njit
def interp_r2(j, x, u, xs, h):
    return u[j-1] + (u[j] - u[j-1]) * (x - xs[j-1]) / h

# @njit
def SI_r2(j, u):
    return (u[j] - u[j-1])**2

# @njit
def CLF(u, xs, dx, fact = 0.2):
    # Set dt using velocity of characteristics
    vel = xs * np.sin(u / xs/xs) * np.cos(u / xs/xs)
    v = np.max(np.array([-np.min(vel), np.max(vel)]))

    # CFL condition
    dt = fact * dx / np.abs(v)
    if dt > 0.01*dx:
       dt = 0.01*dx  # largest timestep allowed
    return dt

# @njit
# def CLF_np(u, xs, dx, fact = 0.45):
#     # Set dt using velocity of characteristics
#     vel = xs * np.sin(u / xs**2) * np.cos(u / xs**2)
#     v_abs = np.abs(vel).max()

#     # Compute dt
#     dt = fact * dx / v_abs
#     if dt > 0.02:  # largest timestep allowed
#         dt = 0.02
#     return dt

# dt_np = CLF_np(u_p[nghost:-nghost], x_phys, h)
# dt    = CLF(u_p[nghost:-nghost], x_phys, h)
# np.testing.assert_almost_equal(dt, dt_np)

# @njit
def der(f, h):
    f_x = np.zeros_like(f)
    f_x[0:3]   = (-147*f[0:3]+360*f[1:4]-450*f[2:5]+400*f[3:6]-225*f[4:7]+72*f[5:8]-10*f[6:9])/(60*h)
    f_x[-3:] = (10*f[-9:-6]-72*f[-8:-5]+225*f[-7:-4]-400*f[-6:-3]+450*f[-5:-2]-360*f[-4:-1]+147*f[-3:])/(60*h)
    f_x[3:-3]  = (-f[0:-6]+9*f[1:-5]-45*f[2:-4]+45*f[4:-2]-9*f[5:-1]+f[6:])/(60*h)
    return f_x



T_final = 2 #8.5*mass**2 + 40*mass

r0 = 15
m = 5

xM = 50.0
h = 0.1

r = 2
nghost = 1

# print('r = ', r)
# print('h = ', h)
# print('T_final =', T_final)

t1 = time.time()

# Grid
nx = int(xM / h + 1 + 2*nghost)
xs = np.array([0.0001*h + (i-nghost)*h for i in range(nx)])
x_phys = xs[nghost:-nghost]

# Initial Data
u_p = np.zeros_like(xs)
u_p[nghost:-nghost] = initial_data(x_phys, m, r0, h)
u_p = BC(u_p, nghost)

t = 0
# print("Starting simulation...")
while t < T_final:
    # Compute dt
    dt = CLF(u_p[nghost:-nghost], x_phys, h)

    # Advance time
    if t + dt > T_final:
        dt = T_final - t
        t = T_final
    else:
        t = t + dt

    # RK timestep
    u = TDV_RK(r, u_p, xs, h, dt, nghost)

    # print(f"Completed {np.round(t*100/T_final, 3) :<3}%\r", end = '')
   
    if t < T_final:
        u_p = u.copy()

# print(f'\nDone in {time.time() - t1} s')
# 
# wenor2 = {}
# wenor2['r'] = r
# wenor2['dx'] = h
# wenor2['x']  = x_phys
# wenor2['t']  = t
# wenor2['B']  = u[nghost:-nghost].copy()
# wenor2['rho_dx'] = der(x_phys**3 * np.sin(u[nghost:-nghost] / x_phys**2)**2, h) / (8 * np.pi * x_phys**2)
# wenor2['rho_dt'] = (u_p[nghost:-nghost] - u[nghost:-nghost]) / x_phys**2 / (4*np.pi*dt)