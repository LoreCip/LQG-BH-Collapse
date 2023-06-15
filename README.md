# LQG-WENO-Reconstruction

This code implements the original formulation of the Weighted Essentially Non-Oscillatory [1] (WENO) to solve the Loop Quantum Gravity equations of a collapsing star under the assumpions of isotropy and uniformity of the matter field.
It solves the equation

$$
\begin{cases}
&\partial_t B(x, t) = - \partial_x  \left[ \frac{1}{2} x^3 \sin^2 \left( \frac{B}{x^2} \right) \right] +  \frac{1}{2}\varepsilon^b \\
&\partial_t \varepsilon^b (x, t) = - \partial_x(\varepsilon^b)\sin\left( 2\frac{B}{x^2} \right)
\end{cases}
$$

for the following initial condition

$$
\begin{cases}
&B(x, 0) = - \frac{x^2}{2} \arccos\left( 1 - 4\frac{M(x)}{x^3} - 2 \frac{\varepsilon^b}{x^2} \right), \qquad M(x) = 4\pi\int_0^x \rho(y)y^2 dy  \\
&\varepsilon^b(x,0) = - \frac{x^2}{a_0^2} \theta(x - r_0) - \frac{r_0^2}{a_0^2} \left(1 - \theta(x - r_0)\right)
\end{cases}
$$

where

$$
\rho(x) = \frac{m}{\frac{4}{3}\pi r_0^3} \left( 1 - \theta(x - r_0)  \right), \quad \theta = \begin{cases} 
0 & \text{ if } x \lt r_0\\ 
0.5 & \text{ if } x = r_0\\ 
1 & \text{ if } x \gt r_0 
\end{cases}
$$

and $m = 5$, and fixing the ratio $\frac{r_0}{a_0} = 0.375$.

## How to install

To install the Fortran code run the bash script **make.sh** using the following command
```bash
./make.sh
```
It creates a folder named _build_ that contains a copy of the source code as present at the moment of compilation, the Makefile and the executable named **run**. It will also create a subfolder named _output_ that will be the default target for the output routines and copy the standard parameter file that can be customized. It requires the HDF5 [2] package to be installed and on the PATH environmental variable.

To not use HDF5 compile using the Makefile with the NOHDF5=1 option
```bash
make NOHDF5=1
```
or modify by hand the **make.sh** file. 

The name and location of the _build_ folder can be customized providing the **make.sh** script a path as an argument
```bash
./make.sh path/to/custom/build/folder
```
To remove an existing build folder and all its content the dedicated command is
```bash
./make.sh remove path/to/build/folder
```
To run the code there are two optional arguments to be given:
```bash
./run [path/to/parameter/file.dat] [path/to/output/folder]
```
If they are not given the standard paths will be used. Please note the code assumes the executable is run from inside the build folder; if executing it from elsewhere these two arguments become mandatory.

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

[2] https://www.hdfgroup.org/solutions/hdf5/
