# Drug-Release-Tracking 


## Project Overview

This project focuses on solving a **coupled system of three nonlinear partial differential equations (PDEs)** that describe interacting transport and reaction dynamics. The system models the temporal and spatial evolution of three dependent variables: $u(x,t)$, $v(x,t)$, and $s(x,t)$, governed by diffusion, reaction, and coupling terms.

The governing equations are:

$$
\begin{aligned}
\frac{\partial u}{\partial t} &= D \frac{\partial^2 u}{\partial x^2} + E \frac{\partial^2 s}{\partial x^2} + f(u,v) \\
\frac{\partial v}{\partial t} &= g(u,v) \\
\frac{\partial s}{\partial t} &= \alpha u - \beta s + \gamma \frac{\partial u}{\partial t}
\end{aligned}
$$

where:

* $D$, $E$: diffusion coefficients
* $\alpha$, $\beta$, $\gamma$: coupling constants
* $f(u,v) = -u(u_b - u) + v(v_b - v)$,
* $g(u,v) = u(u_b - u) - v(v_b - v)$: nonlinear reaction terms
* $u_b$, $v_b$: maximum binding capacities


###  Numerical Methods Used

To solve this system, we implemented and compared **four different numerical approaches**:

1. **Finite Element Method (FEM)**

   * Spatial domain discretized using piecewise linear basis functions.
   * Time evolution computed via implicit time-stepping schemes.
   * Offers accuracy and flexibility in handling complex boundary conditions.

2. **Method of Lines (MOL)**

   * Spatial derivatives discretized using finite difference methods.
   * Transforms the PDE system into a large system of ODEs in time.
   * Solved using standard ODE solvers like Runge-Kutta.

3. **Finite Volume Method (FVM)**

   * Based on conservation laws and control volume integration.
   * Ensures flux balance across cell interfaces.
   * Particularly suited for handling discontinuities and conservation properties.

4. **Physics-Informed Neural Networks (PINN)**

   * Deep learning framework that embeds the governing PDEs into the loss function.
   * Learns the solution using neural networks trained on collocation points.
   * Requires no spatial mesh and can generalize well across the domain.


