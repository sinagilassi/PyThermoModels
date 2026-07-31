# Guidance for Interpreting NRTL and UNIQUAC Binary Interaction Parameters

## 1. Objective

When reading binary interaction parameters from a paper, database, simulator, or YAML file, the agent must determine:

1. The activity-coefficient model: **NRTL** or **UNIQUAC**.
2. Whether the reported parameters describe:
   - an interaction-energy difference,
   - a dimensionless NRTL parameter,
   - the logarithm of a UNIQUAC interaction factor, or
   - the UNIQUAC interaction factor itself.
3. The temperature-dependent equation used.
4. The units and sign convention.
5. The ordered component pair \(i\rightarrow j\).

The agent must not interpret coefficients such as \(A_{ij}\), \(B_{ij}\), or \(C_{ij}\) without identifying their corresponding equation.

---

## 2. Fundamental Definitions

| Model | Energy parameter | Dimensionless parameter | Exponential weighting term |
|---|---|---|---|
| NRTL | \(\Delta g_{ij}=g_{ij}-g_{jj}\) | \(\displaystyle \tau_{ij}=\frac{\Delta g_{ij}}{RT}\) | \(\displaystyle G_{ij}=\exp(-\alpha_{ij}\tau_{ij})\) |
| UNIQUAC | \(\Delta u_{ij}=u_{ij}-u_{jj}\) | Not normally defined as \(\Delta u_{ij}/RT\) | \(\displaystyle \tau_{ij}=\exp\left(-\frac{\Delta u_{ij}}{RT}\right)\) |

### Core distinction

| Model | Meaning of \(\tau_{ij}\) |
|---|---|
| NRTL | Dimensionless interaction-energy difference |
| UNIQUAC | Exponential interaction weighting factor |

Therefore:

\[
\tau_{ij}^{\mathrm{NRTL}}
=
\frac{\Delta g_{ij}}{RT}
\]

but:

\[
\tau_{ij}^{\mathrm{UNIQUAC}}
=
\exp\left(-\frac{\Delta u_{ij}}{RT}\right)
\]

The agent must not apply the UNIQUAC exponential definition to NRTL parameters.

---

## 3. Common NRTL Temperature Forms

For NRTL:

\[
\Delta g_{ij}=RT\tau_{ij}
\]

| Form ID | Reported equation | Calculated NRTL parameter | Corresponding energy equation | Expected coefficient units |
|---|---|---|---|---|
| NRTL-1 | \(\Delta g_{ij}=A_{ij}\) | \(\displaystyle \tau_{ij}=\frac{A_{ij}}{RT}\) | \(\Delta g_{ij}=A_{ij}\) | \(A_{ij}\): J/mol |
| NRTL-2 | \(\Delta g_{ij}/R=A_{ij}\) | \(\displaystyle \tau_{ij}=\frac{A_{ij}}{T}\) | \(\Delta g_{ij}=RA_{ij}\) | \(A_{ij}\): K |
| NRTL-3 | \(\tau_{ij}=A_{ij}\) | \(\tau_{ij}=A_{ij}\) | \(\Delta g_{ij}=RTA_{ij}\) | \(A_{ij}\): dimensionless |
| NRTL-4 | \(\displaystyle \tau_{ij}=A_{ij}+\frac{B_{ij}}{T}\) | Use directly | \(\Delta g_{ij}=R(A_{ij}T+B_{ij})\) | \(A_{ij}\): dimensionless; \(B_{ij}\): K |
| NRTL-5 | \(\Delta g_{ij}=A_{ij}+B_{ij}T\) | \(\displaystyle \tau_{ij}=\frac{A_{ij}}{RT}+\frac{B_{ij}}{R}\) | Use directly | \(A_{ij}\): J/mol; \(B_{ij}\): J/(mol·K) |
| NRTL-6 | \(\displaystyle \tau_{ij}=A_{ij}+\frac{B_{ij}}{T}+\frac{C_{ij}}{T^2}\) | Use directly | \(\displaystyle \Delta g_{ij}=R\left(A_{ij}T+B_{ij}+\frac{C_{ij}}{T}\right)\) | \(A_{ij}\): dimensionless; \(B_{ij}\): K; \(C_{ij}\): K² |
| NRTL-7 | \(\displaystyle \tau_{ij}=A_{ij}+\frac{B_{ij}}{T}+C_{ij}\ln T\) | Use directly | \(\Delta g_{ij}=R[A_{ij}T+B_{ij}+C_{ij}T\ln T]\) | \(A_{ij}\): dimensionless; \(B_{ij}\): K; \(C_{ij}\): dimensionless |
| NRTL-8 | \(\displaystyle \tau_{ij}=A_{ij}+\frac{B_{ij}}{T}+C_{ij}\ln T+D_{ij}T\) | Use directly | \(\Delta g_{ij}=R[A_{ij}T+B_{ij}+C_{ij}T\ln T+D_{ij}T^2]\) | Units must make every term in \(\tau_{ij}\) dimensionless |
| NRTL-9 | General polynomial or simulator-specific expression | Evaluate the declared equation | \(\Delta g_{ij}=RT\tau_{ij}\) | Determine from equation metadata |

### NRTL calculation sequence

| Step | Operation |
|---:|---|
| 1 | Evaluate \(\tau_{ij}(T)\) using the declared correlation. |
| 2 | Evaluate \(\tau_{ji}(T)\) separately. |
| 3 | Determine \(\alpha_{ij}\), usually with \(\alpha_{ij}=\alpha_{ji}\). |
| 4 | Calculate \(G_{ij}=\exp(-\alpha_{ij}\tau_{ij})\). |
| 5 | Use \(\tau_{ij}\) and \(G_{ij}\) in the NRTL activity-coefficient equation. |

The exponential belongs to \(G_{ij}\), not to the fundamental definition of NRTL \(\tau_{ij}\).

---

## 4. Common UNIQUAC Temperature Forms

For standard UNIQUAC:

\[
\tau_{ij}
=
\exp\left(-\frac{\Delta u_{ij}}{RT}\right)
\]

and:

\[
\Delta u_{ij}
=
-RT\ln(\tau_{ij})
\]

| Form ID | Reported equation | Calculated UNIQUAC parameter | Corresponding energy equation | Expected coefficient units |
|---|---|---|---|---|
| UNIQUAC-1 | \(\Delta u_{ij}=A_{ij}\) | \(\displaystyle \tau_{ij}=\exp\left(-\frac{A_{ij}}{RT}\right)\) | \(\Delta u_{ij}=A_{ij}\) | \(A_{ij}\): J/mol |
| UNIQUAC-2 | \(\Delta u_{ij}/R=A_{ij}\) | \(\displaystyle \tau_{ij}=\exp\left(-\frac{A_{ij}}{T}\right)\) | \(\Delta u_{ij}=RA_{ij}\) | \(A_{ij}\): K |
| UNIQUAC-3 | \(-\Delta u_{ij}/R=A_{ij}\) | \(\displaystyle \tau_{ij}=\exp\left(\frac{A_{ij}}{T}\right)\) | \(\Delta u_{ij}=-RA_{ij}\) | \(A_{ij}\): K |
| UNIQUAC-4 | \(\ln\tau_{ij}=A_{ij}\) | \(\tau_{ij}=\exp(A_{ij})\) | \(\Delta u_{ij}=-RTA_{ij}\) | \(A_{ij}\): dimensionless |
| UNIQUAC-5 | \(\displaystyle \ln\tau_{ij}=A_{ij}+\frac{B_{ij}}{T}\) | \(\displaystyle \tau_{ij}=\exp\left(A_{ij}+\frac{B_{ij}}{T}\right)\) | \(\Delta u_{ij}=-R(A_{ij}T+B_{ij})\) | \(A_{ij}\): dimensionless; \(B_{ij}\): K |
| UNIQUAC-6 | \(\Delta u_{ij}=A_{ij}+B_{ij}T\) | \(\displaystyle \tau_{ij}=\exp\left[-\frac{A_{ij}+B_{ij}T}{RT}\right]\) | Use directly | \(A_{ij}\): J/mol; \(B_{ij}\): J/(mol·K) |
| UNIQUAC-7 | \(\displaystyle \ln\tau_{ij}=A_{ij}+\frac{B_{ij}}{T}+\frac{C_{ij}}{T^2}\) | Exponentiate the complete expression | \(\displaystyle \Delta u_{ij}=-R\left(A_{ij}T+B_{ij}+\frac{C_{ij}}{T}\right)\) | \(A_{ij}\): dimensionless; \(B_{ij}\): K; \(C_{ij}\): K² |
| UNIQUAC-8 | \(\displaystyle \ln\tau_{ij}=A_{ij}+\frac{B_{ij}}{T}+C_{ij}\ln T\) | Exponentiate the complete expression | \(\Delta u_{ij}=-R[A_{ij}T+B_{ij}+C_{ij}T\ln T]\) | Units must make \(\ln\tau_{ij}\) dimensionless |
| UNIQUAC-9 | General expression for \(\ln\tau_{ij}\) | Exponentiate the complete expression | \(\Delta u_{ij}=-RT\ln\tau_{ij}\) | Determine term-by-term |
| UNIQUAC-10 | General energy correlation \(\Delta u_{ij}(T)\) | \(\displaystyle \tau_{ij}=\exp[-\Delta u_{ij}(T)/(RT)]\) | Use the declared energy equation | Energy units must be consistent with \(R\) |

---

## 5. Side-by-Side Equation Guidance

| Temperature-correlation family | NRTL implementation | UNIQUAC implementation |
|---|---|---|
| Constant energy | \(\displaystyle \tau_{ij}=\frac{A_{ij}}{RT}\) | \(\displaystyle \tau_{ij}=\exp\left(-\frac{A_{ij}}{RT}\right)\) |
| Energy divided by \(R\) | \(\displaystyle \tau_{ij}=\frac{A_{ij}}{T}\) | \(\displaystyle \tau_{ij}=\exp\left(-\frac{A_{ij}}{T}\right)\) |
| Constant dimensionless parameter | \(\tau_{ij}=A_{ij}\) | \(\ln\tau_{ij}=A_{ij}\), followed by \(\tau_{ij}=e^{A_{ij}}\) |
| Two-coefficient correlation | \(\displaystyle \tau_{ij}=A_{ij}+\frac{B_{ij}}{T}\) | \(\displaystyle \ln\tau_{ij}=A_{ij}+\frac{B_{ij}}{T}\) |
| Inverse-temperature polynomial | Direct equation for \(\tau_{ij}\) | Usually an equation for \(\ln\tau_{ij}\) |
| Logarithmic temperature term | Direct term in \(\tau_{ij}\) | Term inside \(\ln\tau_{ij}\), and therefore inside the exponential |
| General correlation | Evaluate \(\tau_{ij}(T)\) directly | Evaluate \(\ln\tau_{ij}(T)\), then exponentiate |

---

## 6. Parameter-Interpretation Rules

| Rule | Required agent behavior |
|---:|---|
| 1 | Never identify an equation solely from coefficient names such as \(A\), \(B\), \(C\), or \(D\). |
| 2 | Confirm whether \(T\) is in kelvin. |
| 3 | Confirm whether energy parameters are in J/mol, cal/mol, kJ/mol, or another unit. |
| 4 | Use an \(R\) value consistent with the energy units. |
| 5 | Preserve component order because \(ij\) and \(ji\) parameters are generally different. |
| 6 | Do not assume \(A_{ij}=A_{ji}\), \(B_{ij}=B_{ji}\), or \(\tau_{ij}=\tau_{ji}\). |
| 7 | For NRTL, distinguish \(\tau_{ij}\) from \(G_{ij}\). |
| 8 | For UNIQUAC, distinguish \(\tau_{ij}\) from \(\ln\tau_{ij}\). |
| 9 | Inspect the sign convention before applying an exponential. |
| 10 | Reject or flag parameter sets when the equation form or units cannot be established. |

---

## 7. Validation Checks

| Check | NRTL | UNIQUAC |
|---|---|---|
| Dimensionless result | \(\tau_{ij}\) and \(G_{ij}\) must be dimensionless | \(\ln\tau_{ij}\) and \(\tau_{ij}\) must be dimensionless |
| Positivity | \(G_{ij}>0\) | \(\tau_{ij}>0\) |
| Exponential argument | \(-\alpha_{ij}\tau_{ij}\) must be dimensionless | \(-\Delta u_{ij}/RT\) must be dimensionless |
| Pair direction | Calculate \(ij\) and \(ji\) separately | Calculate \(ij\) and \(ji\) separately |
| Temperature range | Must remain within the regression range | Must remain within the regression range |
| Numerical plausibility | Flag overflow or extreme \(G_{ij}\) values | Flag overflow, underflow, or extreme \(\tau_{ij}\) values |

A negative UNIQUAC \(\tau_{ij}\) is always invalid because an exponential is strictly positive.

A negative NRTL \(\tau_{ij}\) can be valid because it is a dimensionless energy difference.
