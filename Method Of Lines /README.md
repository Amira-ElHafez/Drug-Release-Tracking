# Drug Release Tracking Model - Method of Lines (MOL) Implementation

This README provides an overview of a Python-based drug release tracking model implemented using the **Method of Lines (MOL)** to solve a system of three partial differential equations (PDEs). The model simulates the dynamics of unbound drug concentration \( u(x,t) \), bound drug concentration \( v(x,t) \), and polymer stress \( \sigma(x,t) \) in a one-dimensional spatial domain.

---

## 1. Method of Lines (MOL) Overview

The Method of Lines is a numerical technique for solving PDEs by discretizing the spatial domain into a finite grid, transforming the PDEs into a system of ordinary differential equations (ODEs) in time. In this model, MOL is applied using a fourth-order finite difference approximation to discretize spatial derivatives, and the resulting ODE system is solved using SciPy's `odeint` solver.

The model describes:

- **Unbound drug \( u(x,t) \)**: Governed by a diffusion-reaction PDE with coupling to stress.
- **Bound drug \( v(x,t) \)**: Modeled by a reaction ODE with no spatial derivatives.
- **Polymer stress \( \sigma(x,t) \)**: Governed by a PDE with diffusion and coupling to the unbound drug.

---

## 2. Model Parameters, Grid, and Spatial Discretization

This section outlines the physical parameters, spatial grid setup, and the finite difference method used for spatial discretization.

### A) Model Parameters

```python
# Model Parameters
D = 0.6         # Diffusion coefficient for unbound drug u
E = 0.2         # Coupling coefficient with stress σ
ALPHA = 0.2     # Source term coefficient in σ-equation
BETA = 1.0      # Decay coefficient in σ-equation
GAMMA = 1.0     # Coupling coefficient with du/dt in σ-equation
KR = 1.0        # Boundary reaction rate
UA = 0.0        # Ambient drug concentration
UB = 1.0        # Baseline u concentration
VB = 1.0        # Baseline v concentration
```

**Explanation:**

These parameters define the physical behavior of the system:

- \( D \): Controls diffusion of the unbound drug.
- \( E \): Couples the stress field to the unbound drug dynamics.
- \( ALPHA, BETA, GAMMA \): Govern the stress dynamics and its interaction with the drug.
- \( KR, UA \): Define Robin boundary conditions for \( u \).
- \( UB, VB \): Baseline concentrations for reaction terms.

### B) Spatial and Temporal Grid

```python
# Spatial Grid
NX = 26
XL = -0.5
XU = 0.5
xg = np.linspace(XL, XU, NX)

# Temporal Grid
T0 = 0
TF = 2
NOUT = 6
tout = np.linspace(T0, TF, NOUT)
```

**Explanation:**

- The spatial domain \( x \in [-0.5, 0.5] \) is discretized into `NX = 26` uniformly spaced nodes, with grid spacing \( dx = (XU - XL) / (NX - 1) \).
- The temporal domain \( t \in [0, 2] \) is evaluated at `NOUT = 6` output points for visualization, though `odeint` uses adaptive time stepping internally.

### C) Spatial Discretization with Finite Differences

```python
def dss004(xl, xu, n, u):
    dx = (xu - xl) / (n - 1)
    ux = np.zeros(n)
    # Interior points (4th order central difference)
    for i in range(2, n-2):
        ux[i] = (-u[i+2] + 8*u[i+1] - 8*u[i-1] + u[i-2]) / (12*dx)
    # Boundary points (forward/backward differences)
    ux[0] = (-25*u[0] + 48*u[1] - 36*u[2] + 16*u[3] - 3*u[4]) / (12*dx)
    ux[1] = (-3*u[0] - 10*u[1] + 18*u[2] - 6*u[3] + u[4]) / (12*dx)
    ux[n-2] = (-u[n-5] + 6*u[n-4] - 18*u[n-3] + 10*u[n-2] + 3*u[n-1]) / (12*dx)
    ux[n-1] = (3*u[n-5] - 16*u[n-4] + 36*u[n-3] - 48*u[n-2] + 25*u[n-1]) / (12*dx)
    return ux
```

**Explanation:**

The `dss004` function computes the first spatial derivative using a **fourth-order finite difference scheme**:

- **Interior points**: Use a five-point central difference stencil for high accuracy.
- **Boundary points**: Use one-sided fourth-order approximations to maintain consistency.
- Second derivatives (e.g., \( u*{xx}, \sigma*{xx} \)) are computed by applying `dss004` twice.

This transforms the PDEs into a system of \( 3 \times NX \) ODEs for \( u \), \( v \), and \( \sigma \).

---

## 3. Boundary Conditions, Initial Conditions, and Reaction Terms

### A) Robin Boundary Conditions for \( u(x,t) \)

```python
# Apply boundary conditions in the ODE system
ux[0] = -(KR/D) * (UA - u[0])
ux[NX-1] = (KR/D) * (UA - u[NX-1])
sx[0] = 0.0
sx[NX-1] = 0.0
```

**Explanation:**

- **Unbound drug \( u \)**: Robin boundary conditions are applied:
  \[
  D \frac{\partial u}{\partial x} = -KR (UA - u) \text{ at } x = -0.5, \quad D \frac{\partial u}{\partial x} = KR (UA - u) \text{ at } x = 0.5
  \]
  These are incorporated by modifying the first derivative at the boundaries.
- **Stress \( \sigma \)**: Zero Neumann conditions (\( \frac{\partial \sigma}{\partial x} = 0 \)) are applied at both boundaries.

### B) Initial Conditions

```python
def setup_initial_conditions():
    u0 = np.zeros(3 * NX)
    for ix in range(NX):
        u0[ix] = 0.75           # Initial u concentration
        u0[ix + NX] = 0.25      # Initial v concentration
        u0[ix + 2*NX] = 0.0     # Initial s concentration
    return u0
```

**Explanation:**

At \( t = 0 \):

- \( u(x,0) = 0.75 \): Uniform initial unbound drug concentration.
- \( v(x,0) = 0.25 \): Uniform initial bound drug concentration.
- \( \sigma(x,0) = 0 \): Zero initial polymer stress.
  The state vector is concatenated as \( [u, v, \sigma] \), with length \( 3 \times NX \).

### C) Reaction Terms

```python
def f_u(ub, vb, u, v):
    return 0.1 * (ub - u) * (vb - v)

def g_v(ub, vb, u, v):
    return 0.05 * (u - ub) + 0.02 * (v - vb)
```

**Explanation:**

- \( f_u(u, v) \): Source/sink term for the unbound drug PDE, modeling binding/unbinding kinetics.
- \( g_v(u, v) \): Reaction term for the bound drug ODE, representing the rate of binding.
  These terms couple \( u \) and \( v \), with \( UB, VB \) as baseline concentrations.

---

## 4. ODE System and Time Integration

### A) ODE System Definition

```python
def drug_1(U, t):
    u = U[:NX]
    v = U[NX:2*NX]
    s = U[2*NX:3*NX]

    ux = dss004(XL, XU, NX, u)
    sx = dss004(XL, XU, NX, s)

    ux[0] = -(KR/D) * (UA - u[0])
    ux[NX-1] = (KR/D) * (UA - u[NX-1])
    sx[0] = 0.0
    sx[NX-1] = 0.0

    uxx = dss004(XL, XU, NX, ux)
    sxx = dss004(XL, XU, NX, sx)

    ut = np.zeros(NX)
    vt = np.zeros(NX)
    st = np.zeros(NX)

    for i in range(NX):
        ut[i] = D * uxx[i] + E * sxx[i] + f_u(UB, VB, u[i], v[i])
        vt[i] = g_v(UB, VB, u[i], v[i])
        st[i] = ALPHA * u[i] - BETA * s[i] + GAMMA * ut[i]

    return np.concatenate([ut, vt, st])
```

**Explanation:**

The `drug_1` function defines the ODE system:

- Splits the state vector \( U \) into \( u, v, \sigma \).
- Computes first and second spatial derivatives using `dss004`.
- Applies boundary conditions.
- Calculates time derivatives:
  - \( \frac{\partial u}{\partial t} = D u*{xx} + E \sigma*{xx} + f_u(u, v) \)
  - \( \frac{\partial v}{\partial t} = g_v(u, v) \)
  - \( \frac{\partial \sigma}{\partial t} = ALPHA u - BETA \sigma + GAMMA \frac{\partial u}{\partial t} \)
- Returns the concatenated derivative vector.

### B) Time Integration with `odeint`

```python
out = odeint(drug_1, u0, tout)
u_xplot, v_xplot, s_xplot = extract_solutions(out)
```

**Explanation:**

- The ODE system is solved using SciPy's `odeint`, which employs an implicit method (LSODA) suitable for stiff systems.
- The solution is evaluated at `tout` time points and split into \( u, v, \sigma \) for analysis and visualization.

---

## 5. Simulation Results and Visualization

### A) Spatial Profiles at Selected Times

```python
def create_plots(u_xplot, v_xplot, s_xplot, ip_mode):
    if ip_mode == 1:
    #     fig, axes = plt.subplots(3, 1, figsize=(10, 12))
    #     for it in range(NOUT):
    #         axes[0].plot(xg, u_xplot[:, it], label=f't={tout[it]:.1f}')
    #         axes[1].plot(xg, v_xplot[:, it], label=f't={tout[it]:.1f}')
    #         axes[2].plot(xg, s_xplot[:, it], label=f't={tout[it]:.1f}')
    #     axes[0].set_xlabel('x')
    #     axes[0].set_ylabel('u(x,t)')
    #     axes[0].set_title('u(x,t) vs x for different times')
    #     axes[0].legend()
    #     axes[1].set_xlabel('x')
    #     axes[1].set_ylabel('v(x,t)')
    #     axes[1].set_title('v(x,t) vs x for different times')
    #     axes[1].legend()
    #     axes[2].set_xlabel('x')
    #     axes[2].set_ylabel('s(x,t)')
    #     axes[2].set_title('s(x,t) vs x for different times')
    #     axes[2].legend()
    #     plt.tight_layout()
    #     plt.show()
```

**Explanation:**

When `IP = 1`, the code generates three plots showing the spatial profiles of \( u(x,t), v(x,t), \sigma(x,t) \) at different time snapshots (\( t = 0, 0.4, 0.8, 1.2, 1.6, 2.0 \)).

### B) Temporal Evolution at Center Point

```python
def create_plots(u_xplot, v_xplot, s_xplot, ip_mode):
    if ip_mode == 2:
........



**Explanation:**

When `IP = 2`, the code plots the temporal evolution of \( u, v, \sigma \) at the center point (\( x \approx 0 \)).

### C) Conservation Analysis

```python
def analyze_conservation(u_xplot, v_xplot, s_xplot):
    dx = (XU - XL) / (NX - 1)
    total_u = np.trapz(u_xplot, dx=dx, axis=0)
    total_v = np.trapz(v_xplot, dx=dx, axis=0)
    total_s = np.trapz(s_xplot, dx=dx, axis=0)
    total_mass = total_u + total_v + total_s
    print(f"Initial total mass: {total_mass[0]:.6f}")
    print(f"Final total mass: {total_mass[-1]:.6f}")
    print(f"Mass change: {((total_mass[-1]-total_mass[0])/total_mass[0]*100):.4f}%")
```

**Explanation:**

The `analyze_conservation` function computes the total integrated concentrations of \( u, v, \sigma \) over the spatial domain at each time point using the trapezoidal rule. It reports the initial and final total mass and the percentage change, helping to verify conservation properties.

---

## 6. Results


### **Figure 1 – Unbound Drug Concentration \( u(x,t) \)**

<p align="center">
  <img src="Simulation Results/u.png" alt="Figure 1 – Unbound Drug" width="650">
</p>

- **Unbound drug \( u(x,t) \)** decreases over time due to diffusion and binding.

---

### **Figure 2 – Bound Drug Concentration \( v(x,t) \)**

<p align="center">
  <img src="Simulation Results/v.png" alt="Figure 2 – Bound Drug" width="650">
</p>

- **Bound drug \( v(x,t) \)** increases initially then plateaus, indicating saturation or equilibrium.

---

### **Figure 3 – Polymer Stress \( \sigma(x,t) \)**

<p align="center">
  <img src="Simulation Results/s.png" alt="Figure 3 – Polymer Stress" width="650">
</p>


- **Stress \( σ(x,t) \)** initially zero, then increasing due to the drug presence.
