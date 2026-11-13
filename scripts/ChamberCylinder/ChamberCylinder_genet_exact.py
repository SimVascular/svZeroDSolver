# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""svZeroDSolver configuration of the Genet et al. (2023) baseline case.

Builds the ChamberCylinder ventricle with the genet23 Table 1 parameters inside
its circulation, using only existing svZeroDSolver blocks:

  PRESSURE(P_at) -> upstream -[mitral PiecewiseValve]-> ventricle (ChamberCylinder)
    -[aortic PiecewiseValve]-> aortic_root (C_ar) -> arterial (R_p, C_d)
    -> distal (R_d) -> PRESSURE(P_vs)

Valves (genet23 Eq 35, from Caruel et al. 2013 Eq 12): piecewise-linear
conductances K_at (filling), K_iso (isovolumic leak) and K_ar (ejection),
realized by PiecewiseValve with R = 1/K. The aortic side has no closed leak in
Eq 36, so its closed resistance is effectively infinite.

Windkessel (genet23 Eq 36b-c), the C-R-C-R ladder
  C_ar Pdot_ar + (P_ar - P_d)/R_p = Q_ao                (aortic-root node)
  C_d  Pdot_d  + (P_d - P_ar)/R_p = (P_vs - P_d)/R_d    (distal node)
is built from BloodVessel blocks, whose capacitor sits at the outlet node:
BloodVessel(R=0, C=C_ar), BloodVessel(R_p, C_d), BloodVessel(R_d, 0).
BloodVesselRC cannot be used here: it is a flow-through element whose capacitor
stores no volume. The upstream vessel (R = C = 0) only exists because valves
connect two blocks.

C_valve (Table 1, 9e-9 m^3/Pa) is not modeled. It is 43x C_ar, has no node in
Caruel's valve/Windkessel topology, and applied at the ventricle it collapses the
isovolumic phase (see ChamberCylinder_provenance.tex, Discrepancies).

Inputs that genet23 does not tabulate (activation nu(t), atrial pressure P_at(t),
Frank-Starling curve n0(e_c)) follow the authors' input-function code
(genet23 Refs 33/6); P_vs is fitted to the Fig 5 end-systolic volume.
"""
import math

# --- genet23 Table 1: valves and Windkessel (SI) ---
K_at_inv = 1.1111e5   # [Pa s/m^3]
K_ar_inv = 7.6923e4   # [Pa s/m^3]
K_iso_inv = 2.0e9     # [Pa s/m^3]
Cvalve = 9.0e-9       # [m^3/Pa] tabulated but not modeled (see docstring)
Rp = 1.35e7           # [Pa s/m^3]
Rd = 1.0949e8         # [Pa s/m^3]
Car = 2.115e-10       # [m^3/Pa]
Cd = 2.0158e-8        # [m^3/Pa]

# --- genet23 Table 1: geometry and material ---
# Table 1 defines the geometry by the chamber volume 65 mL, wall volume 130 mL
# and aspect ratio L/2Ri = 1.5; its radii and length (19.0, 33.0, 57.1 mm) are
# these values rounded (the rounded ones give 64.76 mL and 130.6 mL), so the
# geometry is derived from the defining volumes.
V_INT, V_WALL, ASPECT = 65e-6, 130e-6, 1.5
_RI = (V_INT / (2.0 * math.pi * ASPECT)) ** (1.0 / 3.0)  # 19.03 mm
_L = 2.0 * ASPECT * _RI                                    # 57.10 mm
_RE = math.sqrt(_RI ** 2 + V_WALL / (math.pi * _L))        # 32.97 mm
GEOM = dict(Ri=_RI, Re=_RE, length=_L, alpha_endo=60.0, alpha_epi=-60.0)
MAT = dict(C1=7.0, C2=0.0, C3=700.0, C4=2.0, C5=50.0, C6=4.0, gamma=70.0,
           k_s=1e8, k_0=260e3, mu=70.0)
SIGMA0 = 65e3  # maximum active stress sigma_0 [Pa]
# |e_c dot| enters the BCS decay exactly (Eq 33); the block's default smoothing
# eps = 1e-4 is a Newton aid, reduced here to a numerically inert value.
FV_REG = 1e-10


# --- ChamberCylinder input from flat parameters -------------------------------
# The block takes its wall material and active stress law as "material" and
# "active_stress" objects (shared with ChamberSphere). The scripts keep the
# constants in one flat dict and split it here.
GENET23_KEYS = ("C1", "C2", "C3", "C4", "C5", "C6")
HO_KEYS = {"ho_a": "a", "ho_b": "b", "ho_af": "a4f", "ho_bf": "b4f", "ho_as": "a4s",
           "ho_bs": "b4s", "ho_afs": "afs", "ho_bfs": "bfs", "ho_k": "k"}
BCS_KEYS = {"sigma_max": "sigma_0", "k_s": "k_s", "k_0": "k_0", "mu": "mu",
            "bcs_alpha": "alpha", "fv_reg": "fv_reg"}


def cylinder_vessel(vv, activation=None, formulation=None, tension=None,
                    frank_starling=None):
    """ChamberCylinder vessel entries from flat parameters.

    vv: geometry, gamma, num_elements, density, bulk_modulus and the material
    (C1-C6 for genet23, or ho_* for holzapfel_ogden) and BCS constants
    (sigma_max, k_s, k_0, mu, bcs_alpha, fv_reg). activation: the
    activation_function object; tension: samples {"t", "T", "fourier_modes",
    "ramp"} of a prescribed tension (active stress uniform_unsteady, no
    activation); frank_starling: {"e_c", "n0"}. formulation: volumetric and
    ends. Returns the keys "formulation", "material", "active_stress",
    "zero_d_element_values" and "activation_function" of the vessel.
    """
    vv = dict(vv)
    if any(float(vv.get(k, 0.0)) != 0.0 for k in HO_KEYS):
        material = {"type": "holzapfel_ogden"}
        material.update({HO_KEYS[k]: vv[k] for k in HO_KEYS if k in vv})
    else:
        material = {"type": "genet23"}
        material.update({k: vv[k] for k in GENET23_KEYS if k in vv})
    for k in list(HO_KEYS) + list(GENET23_KEYS):
        vv.pop(k, None)
    if tension is not None:
        active = {"type": "uniform_unsteady", "t": list(tension["t"]), "T": list(tension["T"]),
                  "fourier_modes": tension.get("fourier_modes", 1),
                  "ramp": 1 if tension.get("ramp", False) else 0}
        for k in BCS_KEYS:
            vv.pop(k, None)
    else:
        active = {"type": "genet23"}
        active.update({BCS_KEYS[k]: vv.pop(k) for k in list(BCS_KEYS) if k in vv})
        if frank_starling:
            active["n0_e_c"], active["n0_values"] = list(frank_starling["e_c"]), list(frank_starling["n0"])
    out = {}
    if formulation:
        out["formulation"] = dict(formulation)
    out.update(material=material, active_stress=active, zero_d_element_values=vv)
    if tension is None and activation is not None:
        out["activation_function"] = activation
    return out


def flat_parameters(vessel):
    """Inverse of cylinder_vessel: the flat parameters of a ChamberCylinder
    vessel (zero_d_element_values plus the material and BCS constants)."""
    vv = dict(vessel["zero_d_element_values"])
    mat = vessel["material"]
    if mat["type"] == "holzapfel_ogden":
        vv.update({k: mat.get(v, 0.0) for k, v in HO_KEYS.items()})
    elif mat["type"] == "genet23":
        vv.update({k: mat.get(k, 0.0) for k in GENET23_KEYS})
    else:
        raise ValueError(f"material {mat['type']}: no flat genet23 constants")
    act = vessel["active_stress"]
    if act["type"] == "genet23":
        defaults = dict(sigma_0=0.0, k_s=1e8, k_0=260e3, mu=70.0, alpha=0.0, fv_reg=1e-4)
        vv.update({k: act.get(v, defaults[v]) for k, v in BCS_KEYS.items()})
    return vv

# --- Input functions (authors' code, genet23 Refs 33/6) ---
MMHG = 134.0  # the authors' Pa-per-mmHg conversion

# Activation nu(t): delay + QRS duration + Bazett-corrected QT + final
# repolarization ramp; +-35 1/s plateaus.
ACT_DELAY = 0.10
ACT_QRSD = 0.080
ACT_QTC = 0.330
ACT_RAMP = 0.010
ACT_AMAX, ACT_AMIN = 35.0, -35.0


def activation_times(period=0.8):
    """Return (tsys, tdias): the nu = 0 upstroke (depolarization time = delay +
    ramp) and downstroke (depolarization end) for a heart period [s]. The
    `nagumo` activation function builds the waveform knees from these, the QRS
    width and the final ramp."""
    QT = ACT_QTC * math.sqrt(period)  # Bazett
    ST = QT - ACT_QRSD
    tsys = ACT_RAMP + ACT_DELAY
    tdias = ACT_QRSD + tsys + ST
    return tsys, tdias


# Atrial pressure P_at(t): 7 mmHg peak, 4.5 mmHg floor, atrial kick timed by the
# PQ interval and the AV-node delay.
PAT_MAX = 7.0 * MMHG            # 938 Pa
PAT_MIN = PAT_MAX - 2.5 * MMHG  # 603 Pa
ATR_PQ = 0.14
ATR_AVNODE = 0.080
# Systemic venous pressure P_vs: not tabulated in genet23 (it only appears in
# Eq 36c). Fitted to the Fig 5 end-systolic volume (74 mL); Caruel Table 1 uses
# 1000 Pa. It sets the afterload floor (~1.8 mL ESV per 500 Pa); EDV, peak
# pressure and torsion are insensitive to it.
P_VS = 350.0


def atrial_pressure(period=0.8):
    """Atrial pressure waveform: PAT_MAX from end-diastole into early systole, a
    linear drop to PAT_MIN over [0.05, 0.10] s, the floor through mid-diastole,
    then the kick PAT_MIN -> PAT_MAX over [T - PQ + ramp, T - AV + ramp].
    Returns (times [s], pressures [Pa])."""
    kick_on = period - ATR_PQ + ACT_RAMP       # 0.670 s at T = 0.8 s
    kick_end = period - ATR_AVNODE + ACT_RAMP  # 0.730 s at T = 0.8 s
    t = [0.0, 0.05, 0.10, kick_on, kick_end, period]
    P = [PAT_MAX, PAT_MAX, PAT_MIN, PAT_MIN, PAT_MAX, PAT_MAX]
    return t, P


TSYS_08, TDIAS_08 = activation_times(0.8)
PAT_KICK_T, PAT_KICK_P = atrial_pressure(0.8)


def build(P_at=PAT_MAX, P_vs=P_VS, sigma_max=SIGMA0, bcs_alpha=12.0, ne=12,
          tsys=TSYS_08, tdias=TDIAS_08, alpha_max=ACT_AMAX, alpha_min=ACT_AMIN,
          act_qrs=ACT_QRSD, act_ramp=ACT_RAMP,
          integrator="stiff", rho_infty=0.5, ncycle=8,
          aortic_Rmax=None, atrial_kick=True, density=1000.0):
    """Return the svZeroDSolver configuration (dict) of the genet23 baseline.

    ne: finite elements through the wall; ncycle: cardiac cycles (the last is
    output). Time step 1 ms (genet23 Sec 3.1). The force-velocity term
    (bcs_alpha > 0) needs the L-stable "stiff" integrator; genet23's own
    energy-consistent midpoint scheme is not available in svZeroDSolver.
    """
    aortic_Rmax = aortic_Rmax if aortic_Rmax is not None else 1.0e12
    vv = dict(GEOM)
    vv.update(MAT)
    vv.update(dict(sigma_max=sigma_max, num_elements=ne, bcs_alpha=bcs_alpha,
                   density=density, fv_reg=FV_REG))
    activation = {"type": "nagumo", "tsys": tsys, "tdias": tdias,
                  "qrs": act_qrs, "ramp": act_ramp,
                  "alpha_max": alpha_max, "alpha_min": alpha_min}
    cfg = {
        "boundary_conditions": [
            {"bc_name": "ATRIUM", "bc_type": "PRESSURE",
             "bc_values": ({"P": list(PAT_KICK_P), "t": list(PAT_KICK_T)}
                           if atrial_kick else {"P": [P_at, P_at], "t": [0.0, 0.8]})},
            {"bc_name": "VENOUS", "bc_type": "PRESSURE",
             "bc_values": {"P": [P_vs, P_vs], "t": [0.0, 0.8]}},
        ],
        "simulation_parameters": {
            "number_of_cardiac_cycles": ncycle,
            # 801 points = 800 steps of exactly 1 ms (genet23 Sec 3.1)
            "number_of_time_pts_per_cardiac_cycle": 801,
            "cardiac_period": 0.8, "steady_initial": False,
            "output_variable_based": True, "absolute_tolerance": 1e-9,
            "maximum_nonlinear_iterations": 50, "output_all_cycles": False,
            "integrator": integrator, "rho_infty": rho_infty,
        },
        # Unloaded reference geometry at rest (rho = eps = 0), so the consistent
        # cavity volume is pi Ri^2 L. genet23 preloads statically to P_at and does
        # not state the arterial initial pressures; starting all pressures at
        # P_at, the arterial transient (R_d C_d = 2.2 s) is periodic by cycle ~12.
        "initial_condition": {"pressure_all": P_at,
                              "volume:ventricle": math.pi * vv["Ri"] ** 2 * vv["length"]},
        "vessels": [
            {"boundary_conditions": {"inlet": "ATRIUM"}, "vessel_id": 0,
             "vessel_length": 1.0, "vessel_name": "upstream",
             "zero_d_element_type": "BloodVessel",
             "zero_d_element_values": {"R_poiseuille": 0.0, "C": 0.0}},
            {"boundary_conditions": {}, "vessel_id": 1, "vessel_length": 1.0,
             "vessel_name": "ventricle", "zero_d_element_type": "ChamberCylinder",
             **cylinder_vessel(vv, activation)},
            # Eq 36b: C_ar at the aortic-valve outlet node P_ar.
            {"boundary_conditions": {}, "vessel_id": 2, "vessel_length": 1.0,
             "vessel_name": "aortic_root", "zero_d_element_type": "BloodVessel",
             "zero_d_element_values": {"R_poiseuille": 0.0, "C": Car}},
            # Eq 36b-c: R_p from P_ar to P_d, C_d at P_d.
            {"boundary_conditions": {}, "vessel_id": 3, "vessel_length": 1.0,
             "vessel_name": "arterial", "zero_d_element_type": "BloodVessel",
             "zero_d_element_values": {"R_poiseuille": Rp, "C": Cd}},
            # Eq 36c: R_d from P_d to the venous pressure P_vs.
            {"boundary_conditions": {"outlet": "VENOUS"}, "vessel_id": 4,
             "vessel_length": 1.0, "vessel_name": "distal",
             "zero_d_element_type": "BloodVessel",
             "zero_d_element_values": {"R_poiseuille": Rd, "C": 0.0}},
        ],
        "valves": [
            {"type": "PiecewiseValve", "name": "mitral",
             "params": {"Rmax": K_iso_inv, "Rmin": K_at_inv,
                        "upstream_block": "upstream", "downstream_block": "ventricle"}},
            {"type": "PiecewiseValve", "name": "aortic",
             "params": {"Rmax": aortic_Rmax, "Rmin": K_ar_inv,
                        "upstream_block": "ventricle", "downstream_block": "aortic_root"}},
        ],
        # Junction names must be "J<n>" for the svZeroDVisualization tools.
        "junctions": [
            {"junction_name": "J0", "junction_type": "NORMAL_JUNCTION",
             "inlet_blocks": ["aortic_root"], "outlet_blocks": ["arterial"]},
            {"junction_name": "J1", "junction_type": "NORMAL_JUNCTION",
             "inlet_blocks": ["arterial"], "outlet_blocks": ["distal"]},
        ],
    }
    return cfg
