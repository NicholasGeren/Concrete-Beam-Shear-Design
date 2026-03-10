# Concrete Beam Shear Design App (ACI 318-19)

A web-based structural engineering tool built with **Python, Streamlit, and GitHub** that analyzes shear forces in reinforced concrete beams and performs shear design checks according to **ACI 318-19**.

The app generates **individual shear force diagrams**, a **superimposed shear envelope**, and calculates required **shear reinforcement spacing**.

---

# Features

### Structural Analysis
- Computes support reactions for multiple loading cases
- Generates shear force diagrams for:
  - Case 1 – Dead + Live load over full span
  - Case 2 – Dead load full span + Live load left half
  - Case 3 – Dead load full span + Live load right half
- Produces a **superimposed shear envelope**
- Highlights **maximum and minimum envelope segments**

### ACI 318-19 Shear Design
Implements the one-way shear provisions of **ACI 318-19 Chapter 22**.

Design condition:

\[
\phi V_n \ge V_u
\]

Where

\[
V_n = V_c + V_s
\]

The application calculates:

### Concrete shear strength \(V_c\)

Based on **ACI 318-19 Table 22.5.5.1**

Option (a):

\[
V_c = 2\lambda \sqrt{f'_c} b_w d
\]

Option (b):

\[
V_c = 0.17\lambda \sqrt{f'_c} b_w d
\]

Option (c):

\[
V_c = 0.66\lambda_s \lambda \rho_w^{1/3} \sqrt{f'_c} b_w d
\]

Where

\[
\lambda_s = \frac{2}{1 + d/10}
\]

### Shear reinforcement contribution

\[
V_s = \frac{A_v f_y d}{s}
\]

The app computes the **required stirrup spacing \(s\)** for the selected reinforcement.

### Critical Section

Shear is evaluated at

\[
x = d
\]

from the face of the support per ACI provisions.

---

# Example Problem

The application was initially developed to solve a beam problem similar to **Problem 6-7**.

Given:
Span = 25 ft
Dead load = 1.8 k/ft
Live load = 1.4 k/ft


Factored load:

\[
w_u = 1.2D + 1.6L
\]

\[
w_u = 4.40 \text{ k/ft}
\]

Support reactions:

\[
R_A = R_B = \frac{wL}{2} = 55 \text{ kips}
\]

Shear at critical section:

\[
V(d) = R_A - w_u d
\]

\[
V(d) \approx 48.6 \text{ kips}
\]

---

# Shear Envelope Visualization

The app generates:

### Individual diagrams
- Case 1
- Case 2
- Case 3

### Envelope diagram
Shows:

- Maximum shear envelope
- Minimum shear envelope
- Governing load case segments

---

# Installation

Clone the repository.

```bash
git clone https://github.com/YOUR_USERNAME/Concrete-Beam-Shear-Design.git
