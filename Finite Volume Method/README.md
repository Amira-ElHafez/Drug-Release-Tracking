
#  Drug Release Tracking Model (Finite Volume Method)

This project simulates a **drug transport and release system** using the **Finite Volume Method (FVM)** for spatial discretization and solves the resulting system of ODEs using Python. It models the dynamics of a **drug diffusing through a tissue-like medium**, binding and unbinding, while being influenced by mechanical stress.

---

##  Model Overview

We solve a coupled system of **three time-dependent PDEs**:

$$
\frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial x^2} + E \frac{\partial^2 s}{\partial x^2} + f(u, v)
$$

$$
\frac{\partial v}{\partial t} = g(u, v)
$$

$$
\frac{\partial s}{\partial t} = \alpha u - \beta s + \gamma \frac{\partial u}{\partial t}
$$

Where:

* $u(x, t)$: free drug concentration
* $v(x, t)$: bound drug concentration
* $s(x, t)$: mechanical stress

### Nonlinear reaction terms:

* $f(u, v) = -u(ub - u) + v(vb - v)$
* $g(u, v) = u(ub - u) - v(vb - v)$

### Parameters:

* `D`, `E`: diffusion coefficients for drug and stress
* `α`, `β`, `γ`: mechanical response constants
* `ub`, `vb`: saturation values
* `kr`: reaction rate coefficient

---

##  What’s Inside

This script uses:

* **FVM** for spatial discretization: conservative and well-suited for flux-based systems
* **Central difference** approximations for second derivatives
* **Custom boundary conditions**:

  * Robin-type at both ends for $u$
  * Neumann (zero-gradient) for $s$
* **Efficient solver**: `solve_ivp` with LSODA (`LSODA`)

---

##  What You Get

When you run the code, it generates:

1. **Spatiotemporal plots** of:

   * Free drug $u(x, t)$
   * Bound drug $v(x, t)$
   * Stress $s(x, t)$

2. **Time series plots** at the center point $x = 0$, showing how each variable evolves

3. **Final profiles** at $t = tf$ for all variables

4. **Tabular data** printed for initial and final values

5. **All figures saved** as PNG files for reporting or analysis


## Running the Code

### Install Requirements

```bash
pip install numpy scipy matplotlib
```

### Run the Script

```bash
python fvm_drug_model.py
```

You’ll see console output and figures will be saved automatically.

---

##  Key Assumptions

* **1D Spatial domain**: $[-0.5, 0.5]$, divided into 26 control volumes
* **Initial conditions**:

  * $u(x, 0) = 0.75$
  * $v(x, 0) = 0.25$
  * $s(x, 0) = 0.0$
* **No-flux conditions** on stress at boundaries
* **Moderate nonlinearity** in reactions—no stiff solver needed in this configuration

---

##  Code Highlights

* Modular design: `dydt` function handles all three PDEs simultaneously
* Function evaluations are counted for benchmarking
* Easy to adjust time range, number of grid points, or parameters
* Easy to switch between **full profile plots** and **center-point evolution** by setting:

```python
ip = 1  # Set to 2 for center point only
```

---

## Educational Value

This model is ideal for students and researchers working on:

* Reaction-diffusion systems
* Coupled PDEs in biomedical engineering
* Numerical methods like FVM
* Drug delivery modeling

It demonstrates how mechanical forces can couple with chemical transport processes, a key concept in **tissue engineering**, **drug delivery**, and **biomechanics**.

---

## Simulation Results
### **Figure 1 – Temporal Evolution at Center Point \( x = 0.02 \) (FVM)**

<p align="center">
  <img src="Simulation%20Results/1.png" alt="Figure 1 – Temporal Evolution of u, v, and s at x = 0.02" width="650">
</p>

<p align="center">
  <em>
    Figure 1 shows the time evolution of the free drug \( u \), bound drug \( v \), and mechanical stress \( s \) at a fixed spatial location \( x = 0.02 \), using the Finite Volume Method (FVM). The top-left panel displays a rapid decay in \( u \), indicating diffusion and reaction depletion. The top-right panel shows \( v \) rising initially and then falling, suggesting binding saturation and unbinding. The bottom-left panel presents the evolution of \( s \), with an initial negative dip followed by gradual recovery, capturing the delayed mechanical response to drug concentration. The bottom-right panel combines all three components, illustrating their coupled behavior over time.
  </em>
</p>

### **Figure 2 – Final Spatial Profiles at \( t = t_f \)**

<p align="center">
  <img src="Simulation%20Results/2.png" alt="Figure 2 – Final Spatial Profiles of u, v, and s" width="650">
</p>

<p align="center">
  <em>
    Figure 2 illustrates the spatial distributions of the free drug \( u(x) \), bound drug \( v(x) \), and mechanical stress \( s(x) \) at the final simulation time \( t = t_f \). The free drug \( u \) exhibits a peak near the center and tapers off toward the boundaries, reflecting symmetric diffusion from the core. The bound drug \( v \) follows a similar profile but with lower amplitude, showing the cumulative effect of binding interactions. The stress \( s \) shows a smooth, possibly parabolic, distribution due to its coupling with the drug concentration and the imposed zero-flux boundary condition. Together, the profiles capture the steady-state behavior resulting from the coupled reaction-diffusion-mechanical system.
  </em>
</p>

