
Finite Element Method (FEM):
----------------------------------------
The Finite Element Method is a powerful numerical technique used for solving partial differential equations (PDEs)
by discretizing the spatial domain into smaller finite elements. In this drug release model, FEM is applied to
the semi-discretized form of the governing equations for unbound drug, bound drug, and polymer stress. The global
mass and stiffness matrices are constructed using linear basis functions, enabling accurate modeling of diffusion,
reaction, and mechanical stress interactions. FEM allows us to incorporate complex boundary conditions (e.g.,
Robin conditions) directly into the system matrices and supports stable integration for stiff PDE systems.

## 1. Model Parameters, Grid, and FEM Matrix Assembly

This section describes how we initialize the physical parameters, create the spatial grid, and construct the **mass** and **stiffness** matrices used in the **Finite Element Method (FEM)**.

---

### A) Model Parameters

```python
# ─── Model Parameters ─────────────────────────────
D = 0.6      # Diffusion coefficient for unbound drug u
E = 0.2      # Coupling coefficient with stress σ
alpha = 0.2  # Source term coefficient in σ-equation
beta = 1.0   # Damping term in σ-equation
gamma = 1.0  # Coupling with du/dt in σ-equation
kr = 1.0     # Robin boundary exchange rate
ua = 0.0     # Boundary drug concentration
```

**Explanation:**

These constants define the physical behavior of the drug transport and polymer stress system. The model incorporates:
- Diffusion and stress coupling via \( D \) and \( E \)
- Dynamic stress response through \( alpha, beta, gamma \)
- Robin-type boundary exchange controlled by \( k_r \)

---

### B) Spatial Grid Definition

```python
# ─── Spatial Grid (Nx = 26, x ∈ [-0.5, 0.5]) ───────
Nx = 26
x_min, x_max = -0.5, 0.5
x = np.linspace(x_min, x_max, Nx)
dx = x[1] - x[0]
```

**Explanation:**

We discretize the spatial domain \( x in [-0.5, 0.5] \) into `Nx = 26` **uniformly spaced nodes**. This creates a structured grid for FEM computation.

- `x`: node coordinates
- `dx`: element width (used to scale FEM basis integrals)

---

### C) FEM Mass and Stiffness Matrix Assembly

```python
# ─── FEM Mass & Stiffness Matrices ────────────────
def fem_matrices(Nx, dx):
    M = np.zeros((Nx, Nx))
    K = np.zeros((Nx, Nx))
    for i in range(Nx - 1):
        Mloc = dx / 6 * np.array([[2, 1], [1, 2]])
        Kloc = 1 / dx * np.array([[1, -1], [-1, 1]])
        for a in range(2):
            for b in range(2):
                M[i+a, i+b] += Mloc[a, b]
                K[i+a, i+b] += Kloc[a, b]
    return M, K
```

**Explanation:**

This function assembles the **global FEM matrices** for linear (1D) basis functions:

Mass Matrix (M): Deals with the time-dependent behavior and storage/inertia. Multiplies the time derivative vector.

Stiffness Matrix (K): Deals with the spatial coupling and resistances due to spatial gradients. Multiplies the vector of nodal values.

### Local Mass Matrix

The local mass matrix \( M^{(e)} \) is:

    M^{(e)} = (dx / 6) * [ [2, 1],
                           [1, 2] ]

### Local Stiffness Matrix

The local stiffness matrix \( K^{(e)} \) is:

    K^{(e)} = (1 / dx) * [ [1, -1],
                           [-1,  1] ]
These local matrices are assembled into global matrices `M` and `K` by looping over each element and adding contributions based on connectivity (adjacent nodes `i+a`, `i+b`).

After construction:

```python
M, K_u = fem_matrices(Nx, dx)
_, K_s = fem_matrices(Nx, dx)
```

- `M`: Global mass matrix shared by `u` and `σ`
- `K_u`: Stiffness matrix for unbound drug `u`
- `K_s`: Stiffness matrix for polymer stress `σ`

---

These matrices are then used in the semi-discrete system of PDEs, where the time derivative terms are multiplied by `M` and the second-order spatial derivatives by `K`.


---


## 2. Boundary Conditions, Initial Conditions, and Reaction Terms

### A) Robin Boundary Conditions for \( u(x,t) \)

```python
# ─── Robin BC for u(x,t) ──────────────────────────
K_u[0, 0] += kr / D
K_u[-1, -1] += kr / D
```

**Explanation:**

We apply **Robin boundary conditions** of the form:

    D * (∂u / ∂x) = ± k_r * (u_a - u)


These are incorporated directly into the **FEM stiffness matrix** \( K_u \) by modifying the diagonal entries at the domain boundaries. This is a standard technique in FEM to weakly impose third-type boundary conditions.

- `kr` is the boundary exchange coefficient.
- `D` is the diffusion coefficient.
- We assume \( u_a = 0 \), which is added later in the RHS vector.

---

### B) Initial Conditions

```python
# ─── Initial Conditions ───────────────────────────
def initial_conditions():
    u0 = np.full(Nx, 0.75)
    v0 = np.full(Nx, 0.25)
    s0 = np.zeros(Nx)
    return np.concatenate([u0, v0, s0])
```

**Explanation:**

At time \( t = 0 \), we initialize the system uniformly:

- \( u(x,0) = 0.75 \): unbound drug is distributed uniformly.
- \( v(x,0) = 0.25 \): bound drug is also uniform but at a lower concentration.
- \(σ(x,0) = 0 \): the polymer stress is initially zero.

These arrays are concatenated into one 1D vector that is used as the initial condition for the time integration routine.

---

### C) Reaction Terms

```python
# ─── Reaction Terms ───────────────────────────────
def f(u, v):
    return - (u * (1 - u) - v * (1 - v))

def g(u, v):
    return (u * (1 - u) - v * (1 - v))
```

**Explanation:**

These functions represent the **binding and unbinding kinetics** of the drug.

- \( f(u,v) \): source term in the PDE for **unbound drug** \( u \), negative of the net binding rate.
- \( g(u,v) \): rate of change for **bound drug** \( v \), positive net binding.

The expressions \( u(1 - u) \) and \( v(1 - v) \) are **logistic terms**, often used to model saturation effects. The idea is:
- As \( u \) increases, binding slows down.
- As \( v \) increases, unbinding becomes more significant.

These nonlinear reaction terms drive the coupling between the \( u \) and \( v \) equations.

---

## 3. Time Integration and Right-Hand Side (RHS) Function

This section implements the **semi-discrete PDE system** (after FEM in space) and integrates it over time using `solve_ivp`.

---

### A) Right-Hand Side Function

```python
# ─── RHS Function for solve_ivp ───────────────────
def rhs(t, y):
    u = y[:Nx]
    v = y[Nx:2*Nx]
    s = y[2*Nx:]

    fu = f(u, v)
    gu = g(u, v)

    rhs_u = -D * (K_u @ u) - E * (K_s @ s) + M @ fu
    rhs_u[0]  += kr * ua
    rhs_u[-1] += kr * ua

    du_dt = np.linalg.solve(M, rhs_u)
    dv_dt = gu
    ds_dt = -beta * s + alpha * u + gamma * du_dt

    return np.concatenate([du_dt, dv_dt, ds_dt])
```

**Explanation:**

This function evaluates the time derivative \( dy/dt ) at each step for the system:

- Splits the solution vector `y` into:
  - `u`: unbound drug
  - `v`: bound drug
  - `s`: stress

- Computes the reaction terms `fu = f(u,v)` and `gu = g(u,v)`
- Assembles the **RHS for the u-equation** including:
  - Diffusion term: \( -D * K_u * u)
  - Stress coupling: \( E * K_s * σ)
  - Reaction source: \( +M f(u,v) \)
  - Robin BC corrections at boundaries

- Solves the linear system \( 
    M * (du/dt) = RHS )
- Evolves `v` and `σ` via their ODEs

Returns the full derivative vector for time integration.

---

### B) Time Integration with `solve_ivp`

```python
# ─── Solve at Snapshot Times ──────────────────────
t_snap = [0.0, 0.4, 0.8, 1.2, 1.6, 2.0]
sol_snap = solve_ivp(rhs, [0, 2], initial_conditions(), t_eval=t_snap,
                     method='BDF', rtol=1e-6, atol=1e-8)

u_snap = sol_snap.y[:Nx, :]
v_snap = sol_snap.y[Nx:2*Nx, :]
s_snap = sol_snap.y[2*Nx:, :]
```

**Explanation:**

We solve the full system over \( t in [0, 2] \) using:

- `method='BDF'`: Backward Differentiation Formula, suitable for stiff systems
- `t_eval=t_snap`: Saves the solution at specific snapshot times for plotting
- `initial_conditions()`: Provides the starting state of `u`, `v`, and `σ`

The output `sol_snap.y` is split back into its components for further analysis and visualization.

---

This setup allows us to track the evolution of drug concentration and stress over time and space using a robust implicit time-stepping scheme.

---
## 3. Simulation Results: Drug and Stress Dynamics

This section presents the results of the simulation:  
1. How the drug concentrations and polymer stress evolve over space and time.  
2. What trends we observe in binding, unbinding, and stress propagation.

---

### A) Spatial Profiles at Selected Times

```python
# Plot u(x,t), v(x,t), σ(x,t) over space at each snapshot
fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
...
plt.show()
```

#### **Output:**

Below are the spatial profiles of the unbound drug, bound drug, and polymer stress at selected time snapshots.

---

### **Figure 1 – Unbound Drug Concentration \( u(x,t) \)**

<p align="center">
  <img src="Numerical/Fig%201.png" alt="Figure 1 – Unbound Drug" width="650">
</p>

- **Unbound drug \( u(x,t) \)** decreases over time as it binds and diffuses outward.

---

### **Figure 2 – Bound Drug Concentration \( v(x,t) \)**

<p align="center">
  <img src="Numerical/Fig%202.png" alt="Figure 2 – Bound Drug" width="650">
</p>

- **Bound drug \( v(x,t) \)** increases initially then plateaus, indicating saturation or equilibrium.

---

### **Figure 3 – Polymer Stress \( \sigma(x,t) \)**

<p align="center">
  <img src="Numerical/Fig%203.png" alt="Figure 3 – Polymer Stress" width="650">
</p>

- **Stress \( σ(x,t) \)** forms a **symmetric bell-like distribution**, initially zero, then increasing due to the drug presence.  

---

### B) Temporal Evolution at \( x = 0 \)

```python
# Solve and plot u(0,t), v(0,t), σ(0,t) over time
...
plt.show()
```

#### **Output:**
### **Figure 4 – Time Evolution at Center Point \( x = 0 \)**

<p align="center">
  <img src="Numerical/Fig%204.png" alt="Figure 4 – Time Evolution at x=0" width="650">
</p>

#### **Interpretation:**

- **Unbound drug \( u(0,t) \)** declines steadily as it diffuses and binds.
- **Bound drug \( v(0,t) \)** rises quickly, then stabilizes.
- **Stress \( σ(0,t) \)** starts from zero, builds up with drug activity, and eventually levels off as the system balances.


---
