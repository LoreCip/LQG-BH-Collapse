import numpy as np

import matplotlib.pyplot as plt

from scipy.integrate import cumulative_trapezoid
from numba import njit, vectorize

@njit
def f(u, x):
    return 0.5 * x*x*x * np.sin(u / x/x)**2 

@njit
def fp(u, x):
    return 0.5 * x * np.sin(2 * u /x/x)

def initial_data(xs, m, r0, h):

    # rho = 3 * m * (1 - heaviside(xs-r0)) / (4 * np.pi * r0**3)
    
    # rho = 3 * m / (4*np.pi * r0**3) * (np.pi/2 - np.arctan(xs - r0))/2
    
    # M = 4*np.pi * cumulative_trapezoid(rho * xs*xs, xs, initial=0)
    # print(f'Mass computed at infinity = {np.round(M[-1],3)}')

    M = m * xs**3 / r0**3 * heaviside(r0 - xs) + m * heaviside(xs - r0) 
        
    B0 = - 0.5 * xs*xs * np.arccos(1 - 4 * M / xs/xs/xs)
    # B0[0] = 0
    return B0

@vectorize
def heaviside(x):
    if x < 0:
        return 0
    else:
        return 1

@njit
def TDV_RK(r, u_prev, fup, xs, h, dt, nghost):
    if r == 2:
        return TDV_RK3(u_prev, xs, fup, h, dt, nghost)
    elif r == 3:
        return TDV_RK4(u_prev, xs, fup, h, dt, nghost)
    else:
        raise Exception("Order not implemented! r = 2 or 3")

@njit
def L(r, j, u, fup, xs, h, nghost):
    """
    Equation 2.7b

    u  --> array of sol
    xs --> array of grid
    
    x_jm = x_{j - 1/2}
    x_jM = x_{j + 1/2}
    """
    x_jM = xs[j] + h / 2
    x_jm = xs[j] - h / 2
    if j == nghost:
        o1 = flux( R(r, j  , x_jM, u, fup, xs, h), R(r, j+1, x_jM, u, fup, xs, h), x_jM) 
        o2 = 0
    else:
        o1 = flux( R(r, j  , x_jM, u, fup, xs, h), R(r, j+1, x_jM, u, fup, xs, h), x_jM) 
        o2 = flux( R(r, j-1, x_jm, u, fup, xs, h), R(r, j  , x_jm, u, fup, xs, h), x_jm)
    return - ( o1 - o2 ) / h

@njit
def R(r, j, x_pt, u, fup, xs, h):
    if r == 2:
        return R2(j, x_pt, u, fup, xs, h)
    elif r == 3:
        return R3(j, x_pt, u, fup, xs, h)

@njit
def alpha(r, j, i, u, fup, x, h):
    if r == 2:
        return alpha_r2(j, i, u, fup, x, h)
    elif r == 3:
        return alpha_r3(j, i, u, fup, x, h)

@njit
def interpolants(r, j, x, u, xs, h):
    if r == 2:
        return interp_r2(j, x, u, xs, h)
    elif r == 3:
        return interp_r3(j, x, u, xs, h)

@njit
def SI(r, j, u):
    if r == 2:
        return SI_r2(j, u)
    elif r == 3:
        return SI_r3(j, u)

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
#     u_dec = np.linspace(np.min(a,b), np.max(a,b), 100)
#     if a <= b:
#         return np.min(f(u_dec, x))
#     elif a > b:
#         return np.max(f(u_dec, x))

@njit
def flux(a, b, x):
    ul = a/x/x
    ur = b/x/x
    FL = 0.5 * x**3 * np.sin(ul)**2 
    FR = 0.5 * x**3 * np.sin(ur)**2 
    ll = np.array([FL, FR])
    if ul <= ur:
        return np.min(ll)
    elif ul > ur:
        if ((ur > -np.pi/2) or (ul < -np.pi/2)):
            return np.max(ll)
        else:
            return 0.5 * x**3

@njit
def BC(arr, nghost):
    arr[:nghost] = arr[nghost]
    arr[-nghost:] = arr[-nghost - 1]
    return arr

@njit
def TDV_RK3(u_prev, xs, fup, h, dt, nghost):
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

    # 1**)
    for j in range(len(rr)-1):
        rr[j] = L(2, j+nghost, u_prev, fup, xs, h, nghost)
    u_n1[nghost] = 0
    u_n1[nghost+1:-nghost] = u_prev[nghost+1:-nghost] + dt * rr[1:]
    u_n1 = BC(u_n1, nghost)

    # 2**)
    for j in range(len(rr)-1):
        rr[j] = L(2, j+nghost, u_n1, fup, xs, h, nghost)
    u_n2[nghost] = 0
    u_n2[nghost+1:-nghost] = u_n1[nghost+1:-nghost] + dt * rr[1:]
    u_n2 = BC(u_n2, nghost)

    # 3**)
    u_n12[nghost] = 0
    u_n12[nghost+1:-nghost] = 3 * u_prev[nghost+1:-nghost] / 4 + u_n2[nghost+1:-nghost] / 4
    u_n12 = BC(u_n12, nghost)

    # 4**)
    for j in range(len(rr)-1):
        rr[j] = L(2, j+nghost, u_n12, fup, xs, h, nghost)
    u_n32[nghost] = 0
    u_n32[nghost+1:-nghost] = u_n12[nghost+1:-nghost] + dt * rr[1:]
    u_n32 = BC(u_n32, nghost)

    # 5**)
    u_final[nghost] = 0
    u_final[nghost+1:-nghost] = u_prev[nghost+1:-nghost] / 3 + 2 * u_n32[nghost+1:-nghost] / 3
    return BC(u_final, nghost)
    
@njit
def R2(j, x_pt, u, fup, xs, h):
    """
    Equation 3.16
    """
    aj0 = alpha(2, j, 0, u, fup, x_pt, h)
    aj1 = alpha(2, j, 1, u, fup, x_pt, h)
    d = aj0 + aj1
    w0 = aj0 / d
    w1 = aj1 / d
    out = w0 * interpolants(2, j, x_pt, u, xs, h) + w1 * interpolants(2, j+1, x_pt, u, xs, h)
    return out

@njit
def alpha_r2(j, i, u, fup, x, h, eps = 1e-6):
    """
    Equations 3.17a & 3.17b 
    """
    if fup[j-1] > 0:
        if i == 0:
            return 1 / 2 / (eps + SI(2, j, u))**2
        elif i == 1:
            return 1 / (eps + SI(2, j+1, u))**2
    elif fup[j-1] <= 0:
        if i == 0:
            return 1 / (eps + SI(2, j, u))**2
        elif i == 1:
            return 1 / 2 / (eps + SI(2, j+1, u))**2

@njit
def interp_r2(j, x, u, xs, h):
    return u[j-1] + (u[j] - u[j-1]) * (x - xs[j-1]) / h

@njit
def SI_r2(j, u):
    return (u[j] - u[j-1])**2

@njit
def TDV_RK4(u_prev, xs, fup, h, dt, nghost):
    """
    uj_0 = u_prev[j]
    uj_1 = uj_0 + dt * L(uj_0)
    uj_2 = uj_0 / 2 + u_j1 / 2 - dt * L(u_j1) / 4 + dt * L(u_j2) / 2
    u_j3 = uj_0 / 9 + 2 * u_j1 / 9 + 2 * u_j2 / 3 - dt * L(u_j1) / 9 - dt * L(u_j2) / 3 + dt * L(u_j3)
    u_new = u_j1 / 3 + u_j2 / 3 + u_j3 / 3 + dt * L(u_j2) / 6 + dt * L(u_j3) / 6
    """

    r_j1 = np.zeros(len(u_prev) - 2*nghost)
    r_j2 = np.zeros(len(u_prev) - 2*nghost)
    r_j3 = np.zeros(len(u_prev) - 2*nghost)

    u_j1    = np.zeros_like(u_prev)
    u_j2    = np.zeros_like(u_prev)
    u_j3    = np.zeros_like(u_prev)
    u_final = np.zeros_like(u_prev)

    for j in range(len(r_j1)-2):
        r_j1[j] = L(3, j+nghost, u_prev, fup, xs, h, nghost)
    u_j1[nghost] = 0
    u_j1[nghost+1:-nghost] = u_prev[nghost+1:-nghost] + dt*r_j1[1:]
    u_j1 = BC(u_j1, nghost)

    for j in range(len(r_j2)-2):
        r_j2[j] = L(3, j+nghost, u_j1, fup, xs, h, nghost)
    u_j2[nghost] = 0
    u_j2[nghost+1:-nghost] = u_prev[nghost+1:-nghost] / 2 + u_j1[nghost+1:-nghost] / 2 - dt*r_j1[1:] / 4 + dt * r_j2[1:] / 2
    u_j2 = BC(u_j1, nghost)

    for j in range(len(r_j3)-2):
        r_j3[j] = L(3, j+nghost, u_j2, fup, xs, h, nghost)
    u_j3[nghost] = 0
    u_j3[nghost+1:-nghost] = u_prev[nghost+1:-nghost]/9 + 2*u_j1[nghost+1:-nghost]/9 + 2*u_j2[nghost+1:-nghost]/3 - dt*r_j1[1:]/9 - dt*r_j2[1:]/3 + dt*r_j3[1:]
    u_j3 = BC(u_j1, nghost)

    u_final[nghost] = 0
    u_final[nghost+1:-nghost] = u_j1[nghost+1:-nghost]/3 + u_j2[nghost+1:-nghost]/3 + u_j3[nghost+1:-nghost]/3 + dt*r_j2[1:]/6 + dt*r_j3[1:]/6
    return BC(u_final, nghost)

@njit
def R3(j, x_pt, u, fup, xs, h):
    """
    Equation 3.18
    """
    aj0 = alpha(3, j, 0, u, fup, x_pt, h)
    aj1 = alpha(3, j, 1, u, fup, x_pt, h)
    aj2 = alpha(3, j, 2, u, fup, x_pt, h)
    d = aj0 + aj1 + aj2
    r1 = aj0 / d
    r2 = aj1 / d
    r3 = aj2 / d
    out = r1 * interpolants(3, j, x_pt, u, xs, h) + r2 * interpolants(3, j+1, x_pt, u, xs, h) + r3 * interpolants(3, j+2, x_pt, u, xs, h)
    return out

@njit
def alpha_r3(j, i, u, fup, x, h, eps = 1e-5):
    """
    Equations 3.17a & 3.17b 
    """
    if fup[j-2] > 0:
        if i == 0:
            return 1 / 12 / (eps + SI(3, j, u))**3
        elif i == 1:
            return 1 / 2 / (eps + SI(3, j+1, u))**3
        elif i == 2:
            return 1 / 4 / (eps + SI(3, j+2, u))**3
    elif fup[j-2] <= 0:
        if i == 0:
            return 1 / 4 / (eps + SI(3, j, u))**3
        elif i == 1:
            return 1 / 2 / (eps + SI(3, j+1, u))**3
        elif i == 2:
            return 1 / 12 / (eps + SI(3, j+2, u))**3

@njit
def interp_r3(j, x, u, xs, h):
    p1 = (u[j] - 2*u[j-1] + u[j-2]) * (x - xs[j-1])*(x - xs[j-1]) / 2 / h/h 
    p2 = (u[j] - u[j-2]) * (x - xs[j-1]) / 2 / h
    p3 = u[j-1] - (u[j] - 2*u[j-1] + u[j-2]) / 24

    return p1 + p2 + p3

@njit
def SI_r3(j, u):
    p1 = (u[j-1] - u[j-2])**2
    p2 = (u[j] - u[j-1])**2
    p3 = (u[j] - 2*u[j-1] + u[j-2])**2

    return 0.5*(p1 + p2) + p3

@njit
def CLF(u, xs, dx, fact = 0.45):
    # Set dt using velocity of characteristics
    vel = xs * np.sin(u / xs/xs) * np.cos(u / xs/xs)
    v = np.max(np.array([-np.min(vel), np.max(vel)]))

    # CFL condition
    dt = fact * dx / np.abs(v)
    if dt > 0.01*dx:
       dt = 0.01*dx  # largest timestep allowed
    return dt

@njit
def der(f, h):
    f_x = np.zeros_like(f)
    f_x[0:3]   = (-147*f[0:3]+360*f[1:4]-450*f[2:5]+400*f[3:6]-225*f[4:7]+72*f[5:8]-10*f[6:9])/(60*h)
    f_x[-3:] = (10*f[-9:-6]-72*f[-8:-5]+225*f[-7:-4]-400*f[-6:-3]+450*f[-5:-2]-360*f[-4:-1]+147*f[-3:])/(60*h)
    f_x[3:-3]  = (-f[0:-6]+9*f[1:-5]-45*f[2:-4]+45*f[4:-2]-9*f[5:-1]+f[6:])/(60*h)
    return f_x


from pathlib import Path
import os
Path("./.frames").mkdir(parents=True, exist_ok=True)

T_final = 15 #8.5*mass**2 + 40*mass

r0 = 15

m = 5

xM = 50.0
h = 0.01

r = 2
if r == 2:
    nghost = 1
elif r == 3:
    nghost = 2

# Grid
nx = int(xM / h + 1 + 2*nghost)
xs = np.array([0.01*h + (i-nghost)*h for i in range(nx)])
x_phys = xs[nghost:-nghost]

# Initial Data
u_p = np.zeros_like(xs)
u_p[nghost:-nghost] = initial_data(x_phys, m, r0, h)
u_p = BC(u_p, nghost)

# Derivative of f(u)
fup = fp(u_p[nghost:-nghost], x_phys)
t = 0
counter = 0
it_counter = 0
y_max = -1

print('Starting simulation...')
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
    u = TDV_RK(r, u_p, fup, xs, h, dt, nghost)

    # Compute derivative for next cycle
    fup = fp(u[nghost:-nghost], x_phys)
    

    if (it_counter % 1000 == 0):

        if t < 12:
            rho = der(x_phys**3 * np.sin(u[nghost:-nghost] / x_phys**2)**2, h) / (8 * np.pi * x_phys**2)
        else:
            rho = (u_p[nghost:-nghost] - u[nghost:-nghost]) / x_phys**2 / (4*np.pi*dt)
        
        fig = plt.figure()
        plt.plot(x_phys[5:], rho[5:])
        plt.title(f"Time = {np.round(t, 4)}")

        if y_max < 1.1*np.max(rho[5:]):
            y_max = 1.1*np.max(rho[5:])

        plt.ylim((0, y_max))
        plt.xlim(-0.1, r0*1.5)
        plt.savefig(f'frames/frame{counter:04d}.png')
        plt.close()

        counter += 1
    it_counter += 1

    if t < T_final:
        u_p = u.copy()

    print(f"Time {np.round(t, 4): <7} at iteration {it_counter: <3}...\r", end = '')
print('Done!')

########
## TO DO
########
# 1) check for existance of frame folder and change name
# 2) check for existance of out_video and change name
# 3) input name of out_video
# 4) option to keep frames

print('Creating out_video.mp4')
os.system("ffmpeg -framerate 30 -loglevel quiet -pattern_type glob -i '.frames/*.png' -c:v libx264 -pix_fmt yuv420p out_video.mp4")
os.system("rm -rf ./.frames")