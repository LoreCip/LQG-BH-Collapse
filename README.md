# LQG-WENO-Reconstruction

This code implements the original formulation of the Weighted Essentially Non-Oscillatory [1] (WENO) to solve the Loop Quantum Gravity equations of a collapsing star under the assumpions of isotropy and uniformity of the matter field.
It solves the equation

$$
\partial_t B(x, t) + \partial_x f(B) = 0, \qquad f(B) = \frac{1}{2} x^3 \sin^2 \left( \frac{B}{x^2} \right) 
$$

for the following initial condition

$$
B(x, 0) = - \frac{x^2}{2} \arccos\left( 1 - 4\frac{M(x)}{x^3} \right), \qquad M(x) = 4\pi\int_0^x \rho(y)y^2 dy 
$$

where

$$
\rho(x) = \frac{m}{\frac{4}{3}\pi r_0^3} \left( 1 - \theta(x - r_0)  \right), \quad \theta = \begin{cases} 
0 & \text{ if } x \lt r_0\\ 
0.5 & \text{ if } x = r_0\\ 
1 & \text{ if } x \gt r_0 
\end{cases}
$$

and $m = 5$, $r_0 = 15$.

## How to install

To install the Fortran code run the bash script **make.sh** using the following command
```bash
./make.sh
```
It creates a folder named _build_ that contains a copy of the source code as present at the moment of compilation, the Makefile and the executable named **run**. The name and location of the _build_ folder can be customized providing the **make.sh** script a path as an argument
```bash
./make.sh path/to/custom/build/folder
```

## Description of the physical problem

W.i.P.

## Description of the algorithm

W.I.P.

## Sources

[1] Xu-Dong Liu, Stanley Osher, Tony Chan,
Weighted Essentially Non-oscillatory Schemes,
Journal of Computational Physics,
Volume 115, Issue 1,
1994,
Pages 200-212,
ISSN 0021-9991,
https://doi.org/10.1006/jcph.1994.1187.
