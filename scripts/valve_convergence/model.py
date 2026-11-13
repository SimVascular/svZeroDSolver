"""Physiologic single-ventricle benchmark models for valve robustness tests.

Units are mmHg, ml, s unless the ``unit_scale`` argument is used to rescale
all pressures (and therefore all resistances, compliances and elastances).
``unit_scale=1333.22`` turns the model into CGS (Barye), which is what the
rescaling robustness test uses: the physics is identical, only the numbers
change.
"""

import copy


def lv_afterload(
    period=0.8,
    emax=2.5,
    emin=0.1,
    vrd=10.0,
    vrs=5.0,
    preload=10.0,
    rp=0.05,
    c_art=1.5,
    rd=1.0,
    pd=5.0,
    rmin=0.005,
    rmax=1.0e4,
    steepness=10.0,
    valve_type="ValveTanh",
    valve_params=None,
    n_cycles=8,
    n_pts=200,
    atol=1e-8,
    unit_scale=1.0,
):
    """LA pressure source -> mitral valve -> LV -> aortic valve -> RCR.

    Returns an svZeroDSolver configuration dictionary.
    """
    p = unit_scale  # pressure scale
    # Resistance = pressure / flow, compliance = volume / pressure,
    # elastance = pressure / volume.
    valves = []
    for name, up, down in [
        ("mitral", "PRELOAD", "ventricle"),
        ("aortic", "ventricle", "aorta"),
    ]:
        params = {"upstream_block": up, "downstream_block": down}
        if valve_type == "ValveTanh":
            params.update(
                Rmax=rmax * p, Rmin=rmin * p, Steepness=steepness / p
            )
        elif valve_type == "PiecewiseValve":
            params.update(Rmax=rmax * p, Rmin=rmin * p)
        else:
            params.update(R=rmin * p)
            params.update(valve_params or {})
        valves.append({"type": valve_type, "name": name, "params": params})

    return {
        "simulation_parameters": {
            "number_of_cardiac_cycles": n_cycles,
            "number_of_time_pts_per_cardiac_cycle": n_pts,
            "cardiac_period": period,
            "absolute_tolerance": atol,
            "output_all_cycles": False,
            "output_variable_based": True,
            "steady_initial": False,
            "max_iter_error_to_warning": True,
        },
        "boundary_conditions": [
            {
                "bc_name": "PRELOAD",
                "bc_type": "PRESSURE",
                "bc_values": {"P": [preload * p, preload * p], "t": [0.0, period]},
            },
            {
                "bc_name": "OUTLET",
                "bc_type": "RESISTANCE",
                "bc_values": {"Pd": pd * p, "R": rd * p},
            },
        ],
        "junctions": [],
        "vessels": [
            {
                "boundary_conditions": {"outlet": "OUTLET"},
                "vessel_id": 0,
                "vessel_length": 10.0,
                "vessel_name": "aorta",
                "zero_d_element_type": "BloodVessel",
                "zero_d_element_values": {"C": c_art / p, "R_poiseuille": rp * p},
            }
        ],
        "valves": valves,
        "chambers": [
            {
                "type": "ChamberElastanceInductor",
                "name": "ventricle",
                "values": {
                    "Emax": emax * p,
                    "Emin": emin * p,
                    "Vrd": vrd,
                    "Vrs": vrs,
                    "Impedance": 0.0,
                },
                "activation_function": {
                    "type": "half_cosine",
                    "t_active": 0.25 * period / 0.8,
                    "t_twitch": 0.35 * period / 0.8,
                },
            }
        ],
        "initial_condition": {"Vc:ventricle": 120.0},
    }


def with_valve(cfg, valve_type, params):
    """Replace every valve in ``cfg`` by ``valve_type`` with ``params``."""
    cfg = copy.deepcopy(cfg)
    for v in cfg["valves"]:
        keep = {
            k: v["params"][k]
            for k in ("upstream_block", "downstream_block")
        }
        v["type"] = valve_type
        v["params"] = dict(params, **keep)
    return cfg
