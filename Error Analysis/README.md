## Error Analysis and Numerical Accuracy Comparison
###  Purpose of Error Analysis :
To assess the reliability of different numerical methods used to simulate the drug release model, we perform a relative error comparison between Method of Lines (MOL), **Finite Volume Method (FVM), and **Finite Element Method (FEM) approaches.
Since we do not have a direct analytical solution to the model, MOL is used as the reference method. The goal is to evaluate how closely the FVM and FEM solutions approximate the MOL solution over time.

---
###  Methodology

At each snapshot time \( t in [0, 2] \), we calculate the relative error of each method with respect to the MOL baseline using the formula:

$$
\text{Relative Error} = \left| \frac{\text{Target} - \text{Reference}}{|\text{Reference}| + 10^{-8}} \right| \times 100\%
$$


The average error across the spatial domain is computed at each time point for a comprehensive overview. For the stress term \( s \), the error at \( t = 0 \) is not computed (since it is initially zero in all methods).

The following tables present the relative error percentages of FVM and FEM methods with respect to the reference MOL solution across different time snapshots. This comparison highlights the accuracy of each method in capturing the dynamics of drug concentration (u), chemical interaction (v), and tissue response (σ) over time.

##  Relative Error Comparison Table (vs MOL)

|  Time | Method |   u(%)  |   v(%)  |   σ(%)  |
|-------|--------|---------|---------|---------|
|  0.00 |   FD   |  0.000  |  0.000  |   ---   |
|       |  FVM   |  3.729  |  1.111  |   ---   |
|       |  FEM   |  0.109  |  0.036  |   ---   |
|  0.40 |   FD   |  0.000  |  0.000  |  0.000  |
|       |  FVM   |  3.524  |  1.038  |  0.723  |
|       |  FEM   |  0.106  |  0.041  |  0.077  |
|  0.80 |   FD   |  0.000  |  0.000  |  0.000  |
|       |  FVM   |  3.355  |  0.978  |  0.742  |
|       |  FEM   |  0.104  |  0.046  |  0.080  |
|  1.20 |   FD   |  0.000  |  0.000  |  0.000  |
|       |  FVM   |  3.215  |  0.927  |  0.762  |
|       |  FEM   |  0.102  |  0.051  |  0.083  |
|  1.60 |   FD   |  0.000  |  0.000  |  0.000  |
|       |  FVM   |  3.100  |  0.890  |  0.780  |
|       |  FEM   |  0.100  |  0.055  |  0.086  |
|  2.00 |   FD   |  0.000  |  0.000  |  0.000  |
|       |  FVM   |  3.004  |  0.861  |  0.797  |
|       |  FEM   |  0.098  |  0.059  |  0.088  |

##  Compute Time Table (seconds)

| Method     | Time (s)  |
|------------|-----------|
| FD (MOL)   | 0.3955    |
| FVM        | 0.2308    |
| FEM        | 0.1006    |

---
From the error table:

- MOL (Method of Lines) is treated as the reference, so by definition, it has 0% relative error.
- FEM (Finite Element Method) consistently shows lower relative errors than FVM, especially in u and σ.
- FVM shows higher errors, particularly as time increases.
