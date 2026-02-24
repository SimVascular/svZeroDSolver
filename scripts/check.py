import sympy as sp
import yaml
from pathlib import Path
import re

# Load YAML files
def load_yaml_files(expmat_path, growth_path):
    """Load the two YAML configuration files."""
    with open(expmat_path, 'r') as f:
        expmat_config = yaml.safe_load(f)
    
    with open(growth_path, 'r') as f:
        growth_config = yaml.safe_load(f)
    
    return expmat_config, growth_config


def extract_residuals(config):
    """Extract residuals from config."""
    return config.get('residuals', [])


def extract_helper_functions(config):
    """Extract helper functions as string."""
    return config.get('helper_functions', '')


def parse_sympy_expression(expr_str, symbols_dict):
    """
    Parse a string expression into a SymPy expression.
    
    Args:
        expr_str: String representation of the expression
        symbols_dict: Dictionary of symbol definitions
    
    Returns:
        SymPy expression
    """
    # Replace 'pi' with sp.pi and 'exp' with sp.exp
    expr_str = expr_str.replace('pi', 'sp.pi')
    expr_str = expr_str.replace('exp', 'sp.exp')
    
    # Evaluate in the context of the symbols
    try:
        expr = sp.sympify(expr_str, locals=symbols_dict)
        return expr
    except Exception as e:
        print(f"Error parsing expression: {expr_str}")
        print(f"Error: {e}")
        return None


def check_equation_equivalence(eq_growth, eq_expmat, eq_num, symbols_dict):
    """
    Check if a growth equation reduces to expmat equation when Theta_c=1, Theta_r=1.
    
    Args:
        eq_growth: Growth model equation
        eq_expmat: Expmat model equation
        eq_num: Equation number
        symbols_dict: Dictionary of symbols
    
    Returns:
        Tuple (is_equivalent, details_dict)
    """
    # Substitute Theta_c = 1 and Theta_r = 1
    Theta_c = symbols_dict.get('Theta_c')
    Theta_r = symbols_dict.get('Theta_r')
    
    eq_growth_reduced = eq_growth.subs({Theta_c: 1, Theta_r: 1})
    eq_growth_simplified = sp.simplify(eq_growth_reduced)
    eq_expmat_simplified = sp.simplify(eq_expmat)
    
    # Check equivalence
    difference = sp.simplify(eq_growth_simplified - eq_expmat_simplified)
    is_equivalent = difference == 0
    
    return is_equivalent, {
        'eq_num': eq_num,
        'eq_growth_original': eq_growth,
        'eq_growth_reduced': eq_growth_reduced,
        'eq_growth_simplified': eq_growth_simplified,
        'eq_expmat': eq_expmat_simplified,
        'difference': difference,
        'is_equivalent': is_equivalent
    }


def main():
    # Define paths
    expmat_path = Path('/Users/katrinkoesler/sv0D_forked/scripts/ChamberSphere_expmat.yaml')
    growth_path = Path('/Users/katrinkoesler/sv0D_forked/scripts/ChamberSphere_growth.yaml')
    
    # Load YAML files
    expmat_config, growth_config = load_yaml_files(expmat_path, growth_path)
    
    # Define all symbols
    Theta_c = sp.symbols('Theta_c', real=True, positive=True)
    Theta_r = sp.symbols('Theta_r', real=True, positive=True)
    Pin, Qin, Pout, Qout = sp.symbols('Pin Qin Pout Qout', real=True)
    radius, velo, stress, tau, volume = sp.symbols('radius velo stress tau volume', real=True)
    dPin_dt, dQin_dt, dPout_dt, dQout_dt = sp.symbols('dPin_dt dQin_dt dPout_dt dQout_dt', real=True)
    dradius_dt, dvelo_dt, dstress_dt, dtau_dt, dvolume_dt = sp.symbols('dradius_dt dvelo_dt dstress_dt dtau_dt dvolume_dt', real=True)
    rho, thick0, radius0 = sp.symbols('rho thick0 radius0', real=True, positive=True)
    C0, C1, C2, C3 = sp.symbols('C0 C1 C2 C3', real=True)
    eta, act, act_plus, sigma_max = sp.symbols('eta act act_plus sigma_max', real=True, positive=True)
    
    # Build symbols dictionary
    symbols_dict = {
        'Theta_c': Theta_c, 'Theta_r': Theta_r,
        'Pin': Pin, 'Qin': Qin, 'Pout': Pout, 'Qout': Qout,
        'radius': radius, 'velo': velo, 'stress': stress, 'tau': tau, 'volume': volume,
        'dPin_dt': dPin_dt, 'dQin_dt': dQin_dt, 'dPout_dt': dPout_dt, 'dQout_dt': dQout_dt,
        'dradius_dt': dradius_dt, 'dvelo_dt': dvelo_dt, 'dstress_dt': dstress_dt, 
        'dtau_dt': dtau_dt, 'dvolume_dt': dvolume_dt,
        'rho': rho, 'thick0': thick0, 'radius0': radius0,
        'C0': C0, 'C1': C1, 'C2': C2, 'C3': C3,
        'eta': eta, 'act': act, 'act_plus': act_plus, 'sigma_max': sigma_max,
        'sp': sp, 'pi': sp.pi, 'exp': sp.exp
    }
    
    # Extract residuals
    residuals_expmat = extract_residuals(expmat_config)
    residuals_growth = extract_residuals(growth_config)
    
    print("="*80)
    print("EQUATION EQUIVALENCE CHECK: ChamberSphere_growth vs ChamberSphere_expmat")
    print("="*80)
    print(f"\nWhen Theta_c = 1 and Theta_r = 1")
    print(f"Total equations to check: {len(residuals_expmat)}\n")
    
    # Parse and check each equation
    all_equivalent = True
    results = []
    
    for i, (eq_growth_str, eq_expmat_str) in enumerate(zip(residuals_growth, residuals_expmat)):
        print(f"\n{'='*80}")
        print(f"EQUATION {i + 1}")
        print(f"{'='*80}")
        
        # Parse expressions
        eq_growth = parse_sympy_expression(eq_growth_str, symbols_dict)
        eq_expmat = parse_sympy_expression(eq_expmat_str, symbols_dict)
        
        if eq_growth is None or eq_expmat is None:
            print("SKIPPED: Could not parse equation")
            continue
        
        # Check equivalence
        is_equiv, details = check_equation_equivalence(eq_growth, eq_expmat, i + 1, symbols_dict)
        results.append(details)
        
        print(f"\nGrowth model (original):")
        print(f"  {eq_growth}")
        print(f"\nGrowth model (after substitution Theta_c=1, Theta_r=1):")
        print(f"  {details['eq_growth_reduced']}")
        print(f"\nGrowth model (simplified):")
        print(f"  {details['eq_growth_simplified']}")
        print(f"\nExpmat model (simplified):")
        print(f"  {details['eq_expmat']}")
        print(f"\nDifference (should be 0 if equivalent):")
        print(f"  {details['difference']}")
        print(f"\n✓ EQUIVALENT" if is_equiv else "✗ NOT EQUIVALENT")
        
        all_equivalent = all_equivalent and is_equiv
    
    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    
    equivalent_count = sum(1 for r in results if r['is_equivalent'])
    total_count = len(results)
    
    print(f"\nEquations checked: {total_count}")
    print(f"Equivalent equations: {equivalent_count}")
    print(f"Non-equivalent equations: {total_count - equivalent_count}")
    
    if all_equivalent and total_count > 0:
        print("\n✓ SUCCESS: All equations are mathematically equivalent!")
        print("  ChamberSphere_growth with Theta_c=1 and Theta_r=1 reduces to")
        print("  ChamberSphere_expmat (no-growth model)")
    elif total_count == 0:
        print("\n✗ ERROR: No equations were successfully parsed")
    else:
        print("\n✗ FAILURE: Some equations are NOT equivalent")
        print("\nNon-equivalent equations:")
        for r in results:
            if not r['is_equivalent']:
                print(f"  - Equation {r['eq_num']}")


if __name__ == "__main__":
    main()

