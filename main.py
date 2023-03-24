import numpy as np
from tqdm import tqdm

def main():
    r = 3
    h = 1e-2
    dt = 1e-3    
    N_it = 500

    if r == 2:
        n_ghost = 1
    elif r == 3:
        n_ghost = 2

    #### GRID
    nx = int(2 / h + 1 + 2*n_ghost)
    xs = np.array([(i-n_ghost)*h - 1 for i in range(nx)])

    ### Initial Data
    u_prev = np.ones(len(xs)) # 1 + phys + 1
    # Phys
    u_prev[n_ghost:-n_ghost] = np.sin(4 * np.pi * xs[n_ghost:-n_ghost])
    # Ghosts
    u_prev = PBC(u_prev, n_ghost)

    sol = np.zeros((len(u_prev[n_ghost:-n_ghost]), 2))
    sol[:, 0] = u_prev[n_ghost:-n_ghost].copy()

    t = 0
    for i in tqdm(range(N_it)):
        u_next = TDV_RK(r, u_prev, xs, h, dt, n_ghost)
        u_prev = u_next.copy()

        t = t + dt

    sol[:, 1] = u_next[n_ghost:-n_ghost].copy()


def argparse():
    ## W.i.P.
    return 1

if __name__ == '__main__':
    ## Argparser
    args = argparse()
    ## Main
    main()