#  Drug Release Tracking via Numerical and Deep Learning Methods

📁 **Full Report**: [Click here to view the full PDF report](https://drive.google.com/file/d/1sjhQngD_V0-lfyJuTGpfvhyAA7eJTuBH/view?usp=sharing)
📁 **Slides**: [Click here to view the full PDF Slides]([https://drive.google.com/file/d/1sjhQngD_V0-lfyJuTGpfvhyAA7eJTuBH/view?usp=sharing](https://www.canva.com/design/DAGqlSWfstI/_3BSJV24NPG8SgkdEklldQ/edit?utm_content=DAGqlSWfstI&utm_campaign=designshare&utm_medium=link2&utm_source=sharebutton)

##  Overview

This project investigates and compares four advanced techniques **Method of Lines (MOL)**, **Finite Volume Method (FVM)**, **Finite Element Method (FEM)**, and **Physics-Informed Neural Networks (PINNs)** to solve a nonlinear system of coupled PDEs modeling drug diffusion and mechanical stress in polymeric tissues.

---

## Problem Description

We consider a bio-mathematical model for drug release involving three dependent variables:

- u(x,t) → Unbound (free) drug  
- v(x,t) → Bound drug  
- σ(x,t) → Stress in the polymer matrix  

The model incorporates **nonlinear binding kinetics**, **diffusion**, and **stress coupling**. Governing equations:

$$
\frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial x^2} + E \frac{\partial^2 \sigma}{\partial x^2} + f(u, v)
$$

$$
\frac{\partial v}{\partial t} = g(u, v)
$$

$$
\frac{\partial \sigma}{\partial t} = \alpha u - \beta \sigma + \gamma \frac{\partial u}{\partial t}
$$

---

##  Numerical Methods

###  Method of Lines (MOL)
- Spatial discretization via finite differences  
- Temporal integration using *LSODA*  
- Good accuracy, but relatively *computationally expensive*

###  Finite Volume Method (FVM)
- Volume-integrated discretization with central-difference fluxes  
- Robin boundary conditions implemented  
- Conservative and physically grounded, but slightly less accurate

###  Finite Element Method (FEM)
- Weak form via linear elements  
- Mass and stiffness matrix formulation  
- *Best trade-off* between accuracy and runtime

###  PINNs (Physics-Informed Neural Networks)
- Mesh-free neural network model  
- Trained to minimize residuals of PDEs + initial/boundary loss  
- Promising deep learning alternative, but training-sensitive

---

##  Results Summary

| *Time* | *Method* | *u(%)* | *v(%)* | *σ(%)* |
|---------:|------------|---------:|---------:|---------:|
| 0.00     | MOL        | 0.000    | 0.000    | —        |
|          | FVM        | 3.729    | 1.111    | —        |
|          | FEM        | 0.109    | 0.036    | —        |
| 0.40     | MOL        | 0.000    | 0.000    | 0.000    |
|          | FVM        | 3.524    | 1.038    | 0.723    |
|          | FEM        | 0.106    | 0.041    | 0.077    |

⏱ *Compute Time*:
- MOL: 0.3094 s  
- FVM: 0.1672 s  
- FEM: 0.1151 s  

Note: The plots, results, and performance graphs for each method are provided within their respective method sections in their README. Scroll through to view simulation outputs for MOL, FVM, FEM, PINN and accuracy comparisons.

---

##  PINN Architecture

- 8 hidden layers, 100 neurons each  
- tanh activation, Xavier initialization  
- Adaptive loss weighting schedule  
- Trained with Adam optimizer (lr = 0.0008) for 3000 epochs  

---

##  Conclusions

- *FEM* outperformed other numerical methods in both *accuracy and efficiency*
- *PINNs* are a flexible future path, requiring tuning and more resources  
- FVM shows strong conservation but needs higher-order schemes to match FEM

---

##  Future Directions

- Upgrade MOL with adaptive or higher-order schemes  
- Extend FEM to 2D/3D geometries  
- Introduce TVD/flux-limiter techniques for FVM  
- Improve PINNs via dynamic loss weighting and physics-informed regularization

---

##  References

1. W. E. Schiesser, Differential Equation Analysis in Biomedical Science  
2. M. Raissi et al., Physics-Informed Neural Networks, JCP, 2019  
3. Zhou & Wu, Finite Element Analysis of Drug Diffusion, 2006  
4. LeVeque, Finite Volume Methods for Hyperbolic Problems, 2002  
5. Liu et al., Conservative FVM in Drug Modeling, 2017  
6. J. N. Reddy, Introduction to FEM, 4th ed.  
7. Brenner & Scott, FEM for Coupled PDEs, SIAM, 2007

---
##  Team Members
| Name              | ID           | Contribution                                      |
|-------------------|--------------|---------------------------------------------------|
| Aliaa Mahmoud     |    9230592   | FEM Method, accuracy comparison                   |
| Arwa Mohamed      |  9230203     | ML|
| Amira El-Sayed    |   9230244    | MOL implementation, ML                            |
| Mai Mahmoud       | 9230934      | FEM Method, accuracy comparison                   |
| Mahmoud Abdullah  |  9220795     | ML                                                 |
| Mohamed Mandour   |  9230759     | FVM Method                                      |
| Fady Osama        |  91241293    | FVM Method  , accuracy comparison                |

**Special thanks to Dr. Muhammad Rushdi and TA Alaa Tarek for their guidance.**
