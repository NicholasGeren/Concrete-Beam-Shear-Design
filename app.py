# app.py -- Shear envelope + ACI318-19 shear design
# Run: pip install streamlit numpy matplotlib pandas
#       streamlit run app.py

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

st.set_page_config(page_title="Shear Envelope + ACI 318-19 Shear Design", layout="wide")

st.title("Shear envelope & ACI 318-19 shear design (Problem P6-7 style)")

# -----------------------
# Inputs
# -----------------------
with st.sidebar.form("inputs"):
    st.header("Inputs (defaults from Figure P6-7)")
    span_ft = st.number_input("Span, L (ft)", value=25.0, format="%.3f")
    D_kpf = st.number_input("Dead load D (kips/ft)", value=1.8, format="%.3f")
    L_kpf = st.number_input("Live load L (kips/ft)", value=1.4, format="%.3f")

    # load factors (LRFD style default)
    lf_D = st.number_input("Dead load factor", value=1.2, format="%.3f")
    lf_L = st.number_input("Live load factor", value=1.6, format="%.3f")

    # section
    bw_in = st.number_input("Web width b_w (in)", value=12.0, format="%.3f")
    d_in = st.number_input("Effective depth d (in)", value=17.5, format="%.3f")
    fc = st.number_input("Concrete f'c (psi)", value=4500.0, format="%.1f")

    # reinforcement & stirrups
    rho_w = st.number_input("Longitudinal reinforcement ratio ρ_w (decimal)", value=0.01, format="%.4f",
                             help="As/(b_w * d) if you know As you can set rho_w accordingly")
    fy = st.number_input("Stirrup fy (psi)", value=40000.0, format="%.0f")
    # No.3 double-leg default area = 0.11 in² per leg -> 0.22 in² double-leg
    Av_default = 0.22
    Av = st.number_input("Area of stirrup legs Av (in^2) (per stirrup, total legs)", value=Av_default, format="%.4f",
                         help="Total shear-leg area per stirrup (e.g., #3 double-leg ≈ 0.22 in²)")

    phi_shear = st.number_input("φ for shear", value=0.75, format="%.3f")
    use_aci_table = st.checkbox("Use ACI 318-19 Table 22.5.5.1 (detailed Vc)", value=True)
    rho_w_input_method = st.checkbox("Prefer to provide As (in²) instead of ρ_w", value=False)
    if rho_w_input_method:
        As = st.number_input("As (area of tension longitudinal steel) (in²)", value=1.33, format="%.3f")
        rho_w = As / (bw_in * d_in)

    st.form_submit_button("Update")

# -----------------------
# Helper functions
# -----------------------
def make_loads(span, D_kpf, L_kpf, lf_D, lf_L, case):
    """
    return list of loads in format (w_kips_per_ft, a_ft, b_ft)
    case 1: factored dead+live on entire beam
    case 2: factored dead on entire beam + factored live on left half
    case 3: factored dead on entire beam + factored live on right half
    """
    Df = lf_D * D_kpf
    Lf = lf_L * L_kpf
    loads = []
    if case == 1:
        loads.append((Df + Lf, 0.0, span))
    elif case == 2:
        loads.append((Df, 0.0, span))
        loads.append((Lf, 0.0, span / 2.0))
    elif case == 3:
        loads.append((Df, 0.0, span))
        loads.append((Lf, span / 2.0, span))
    else:
        raise ValueError("unknown case")
    return loads

def reactions_for_loads(span, loads):
    """ Solve RA, RB for simply supported beam under piecewise UDLs """
    total_vert = sum(w * (b - a) for (w, a, b) in loads)
    sum_mom = sum(w * (b - a) * ((a + b) / 2.0) for (w, a, b) in loads)
    RB = sum_mom / span
    RA = total_vert - RB
    return RA, RB

def shear_profile(span, loads, RA, npts=400):
    xs = np.linspace(0, span, npts)
    Vs = []
    for x in xs:
        V = RA
        for (w, a, b) in loads:
            if x <= a:
                continue
            overlap = min(b, x) - a
            overlap = max(0.0, overlap)
            V -= w * overlap
        Vs.append(V)
    return xs, np.array(Vs)

# ACI shear functions
def Vc_approx(fc, bw, d):
    """classic approximate Vc = 2*sqrt(fc')*b_w*d (units -> kips)"""
    return 2.0 * math.sqrt(fc) * bw * d / 1000.0

def lambda_s_for_d(d):
    """ACI 318-19 Eq. 22.5.5.1.3: lambda_s = 2/(1 + d/10) <= 1"""
    val = 2.0 / (1.0 + d / 10.0)
    return min(1.0, val)

def Vc_tablec(fc, bw, d, rho_w, lambda_factor=1.0, Nu_over_Ag=0.0):
    """
    Implement ACI 318-19 Table 22.5.5.1(c) style expression as:
      v_c = 0.33 * lambda_s * lambda * (rho_w)**(1/3) * sqrt(fc)
      Vc = v_c * b_w * d (convert to kips)
    NOTE: ACI table allows different expressions (a,b,c). This is a faithful *implementation
    of the common interpretation* used by design tools and by the ACI updates.
    The app also caps Vc <= 5 * lambda * sqrt(fc) * b_w * d (per 22.5.5.1.1).
    """
    lambda_s = lambda_s_for_d(d)
    v_c = 0.33 * lambda_s * lambda_factor * (rho_w ** (1.0 / 3.0)) * math.sqrt(fc)
    Vc_kips = v_c * bw * d / 1000.0
    Vc_cap = 5.0 * lambda_factor * math.sqrt(fc) * bw * d / 1000.0
    Vc_kips = min(Vc_kips, Vc_cap)
    return Vc_kips, lambda_s, Vc_cap

def required_stirrup_spacing(Vu_kips, Vc_kips, Av_in2, fy_psi, d_in, phi=0.75):
    """Compute required maximum spacing s (in) given Av, fy, d to supply Vs = (Vu/phi - Vc)"""
    V_needed = max(0.0, Vu_kips / phi - Vc_kips)  # kips
    if V_needed == 0:
        return float('inf'), V_needed
    # Vs (lb) = Av * fy * d / s  => s = Av * fy * d / Vs(lb)
    Vs_lb = V_needed * 1000.0
    s_in = Av_in2 * fy_psi * d_in / Vs_lb
    return s_in, V_needed

# -----------------------
# Build cases, diagrams
# -----------------------
cases = [1, 2, 3]
case_profiles = {}
for case in cases:
    loads = make_loads(span_ft, D_kpf, L_kpf, lf_D, lf_L, case)
    RA, RB = reactions_for_loads(span_ft, loads)
    xs, Vs = shear_profile(span_ft, loads, RA, npts=800)
    case_profiles[case] = dict(loads=loads, RA=RA, RB=RB, xs=xs, Vs=Vs)

# build envelope
xs = case_profiles[1]["xs"]  # same x grid for all
Vs_all = np.vstack([case_profiles[c]["Vs"] for c in cases])  # shape (3, n)
Vmax = Vs_all.max(axis=0)
Vmin = Vs_all.min(axis=0)

# find segments where envelope attains maxima/minima
imax_segments = np.where(np.isclose(Vmax, Vs_all.max(axis=0)))[0]  # placeholder
# for coloring, we'll find where each case equals the envelope
envelope_source_max = np.argmax(Vs_all, axis=0)  # 0..2 mapping to case index-1
envelope_source_min = np.argmin(Vs_all, axis=0)

# -----------------------
# Display plots and results
# -----------------------
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("Individual shear diagrams (kips)")
    fig, ax = plt.subplots(figsize=(8, 4))
    colors = {1: "C0", 2: "C1", 3: "C2"}
    for c in cases:
        ax.plot(case_profiles[c]["xs"], case_profiles[c]["Vs"], label=f"Case {c}", color=colors[c])
    ax.axhline(0, color="k", linewidth=0.6)
    ax.set_xlabel("x (ft)")
    ax.set_ylabel("Shear (kips)")
    ax.legend()
    st.pyplot(fig)

with col2:
    st.subheader("Shear envelope (superimposed)")
    fig2, ax2 = plt.subplots(figsize=(8, 4))
    # plot all cases faint
    for c in cases:
        ax2.plot(case_profiles[c]["xs"], case_profiles[c]["Vs"], color="lightgray", linewidth=0.7)
    ax2.plot(xs, Vmax, label="Envelope max", color="tab:red", linewidth=2.0)
    ax2.plot(xs, Vmin, label="Envelope min", color="tab:blue", linewidth=2.0)
    # highlight which case gives envelope (color segments)
    # we color max segments by which case produces them
    for c in cases:
        mask = (envelope_source_max == (c - 1))
        if mask.any():
            ax2.plot(xs[mask], Vmax[mask], linewidth=3.0, label=f"max from case {c}", alpha=0.9)
        mask_min = (envelope_source_min == (c - 1))
        if mask_min.any():
            ax2.plot(xs[mask_min], Vmin[mask_min], linewidth=3.0, linestyle="--", label=f"min from case {c}", alpha=0.9)
    ax2.axhline(0, color="k", linewidth=0.6)
    ax2.set_xlabel("x (ft)")
    ax2.set_ylabel("Shear (kips)")
    ax2.legend(loc="upper right", bbox_to_anchor=(1.35, 1.0))
    st.pyplot(fig2)

# -----------------------
# Shear design at critical section (x = d)
# -----------------------
st.header("Shear design at critical section (x = d from support)")

d_ft = d_in / 12.0
# compute shear V at x=d for each case
V_at_d = {}
for c in cases:
    RA = case_profiles[c]["RA"]
    loads = case_profiles[c]["loads"]
    # compute shear at x = d_ft
    V = RA
    for (w, a, b) in loads:
        if d_ft > a:
            overlap = min(b, d_ft) - a
            overlap = max(0.0, overlap)
            V -= w * overlap
    V_at_d[c] = V

# pick worst Vu (absolute) among locations (left support side)
Vu_kips = max(abs(V_at_d[c]) for c in cases)

st.write(f"Critical section distance: d = {d_in:.2f} in ({d_ft:.4f} ft)")
st.write("Shear at d for each case (kips):")
st.table(pd.DataFrame.from_dict(V_at_d, orient='index', columns=["V_at_d (kips)"]))

st.write(f"Design shear to check (Vu) = max |V(d)| among cases = **{Vu_kips:.3f} kips**")

# Concrete contribution Vc (two options)
Vc_approx_kips = Vc_approx(fc, bw_in, d_in)
Vc_tablec_kips, lambda_s_val, Vc_cap = Vc_tablec(fc, bw_in, d_in, rho_w, lambda_factor=1.0)

st.write("Concrete shear contribution options:")
st.write(f"- Approximate (classic) Vc = 2·√(f'c)·b_w·d  = **{Vc_approx_kips:.3f} kips**")
st.write(f"- ACI 318-19 Table option (c) Vc = **{Vc_tablec_kips:.3f} kips**; size factor λ_s = {lambda_s_val:.3f}; Vc cap = {Vc_cap:.3f} kips")
st.write("Note: if Av ≥ Av_min you may use Eq (a) or (b) from Table 22.5.5.1; else use (c).")

# required spacing for provided Av
s_in, V_needed_kips = required_stirrup_spacing(Vu_kips, Vc_tablec_kips, Av, fy, d_in, phi=phi_shear)
s_in_approx, V_needed_kips_approx = required_stirrup_spacing(Vu_kips, Vc_approx_kips, Av, fy, d_in, phi=phi_shear)

st.write("---")
st.subheader("Stirrup spacing results (No. inputed stirrup Av as area per stirrup)")

st.write(f"Using Table 22.5.5.1(c) Vc: required Vs = {V_needed_kips:.3f} kips -> required spacing s = {s_in:.2f} in")
st.write(f"Using approximate Vc: required Vs = {V_needed_kips_approx:.3f} kips -> required spacing s = {s_in_approx:.2f} in")

# spacing limits (practical)
s_max_code = min(d_in/2.0, 12.0)  # typical ACI spacing limits; you can refine
st.write(f"Typical code spacing limit used (min(d/2, 12 in)) = {s_max_code:.2f} in")

if s_in < 1.0:
    st.warning("Required spacing (based on Table Vc) is extremely tight (<1 in). Re-check reinforcement assumptions or use larger stirrups / additional longitudinal reinforcement.")
if s_in > s_max_code:
    st.info("Provided Av with that spacing would be adequate under the chosen model; verify detailing per ACI maximum spacing limits.")

# Provide raw arrays (optionally downloadable)
st.header("Export: shear arrays")
df_export = pd.DataFrame({
    "x_ft": xs,
    "V_case1_kips": case_profiles[1]["Vs"],
    "V_case2_kips": case_profiles[2]["Vs"],
    "V_case3_kips": case_profiles[3]["Vs"],
    "V_envelope_max_kips": Vmax,
    "V_envelope_min_kips": Vmin
})
st.dataframe(df_export.head(10))
csv = df_export.to_csv(index=False).encode('utf-8')
st.download_button("Download shear arrays CSV", data=csv, file_name="shear_arrays.csv", mime="text/csv")

st.markdown("----")
st.info("References: ACI 318-19 Table 22.5.5.1 (shear Vc expressions and λs size factor); φ for shear = 0.75; stirrup Vs modelled as Av·fy·d/s. See app header for full citations.")
