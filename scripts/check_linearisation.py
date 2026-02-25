import sympy as sp
from sympy import symbols, simplify, pi, exp, Abs, sqrt
import re

# Define all symbols
Theta_c = sp.Symbol('Theta_c', real=True, positive=True)
Theta_r = sp.Symbol('Theta_r', real=True, positive=True)
Pin, Qin, Pout, Qout = symbols('Pin Qin Pout Qout', real=True)
radius, velo, stress, tau, volume = symbols('radius velo stress tau volume', real=True)
dPin_dt, dQin_dt, dPout_dt, dQout_dt = symbols('dPin_dt dQin_dt dPout_dt dQout_dt', real=True)
dradius_dt, dvelo_dt, dstress_dt, dtau_dt, dvolume_dt = symbols('dradius_dt dvelo_dt dstress_dt dtau_dt dvolume_dt', real=True)
rho, thick0, radius0 = symbols('rho thick0 radius0', real=True, positive=True)
C0, C1, C2, C3 = symbols('C0 C1 C2 C3', real=True)
W1, W2 = symbols('W1 W2', real=True, positive=True)
eta, sigma_max = symbols('eta sigma_max', real=True, positive=True)

# Build symbols dictionary for sympify
symbols_dict = {
    'Theta_c': Theta_c, 'Theta_r': Theta_r,
    'Pin': Pin, 'Qin': Qin, 'Pout': Pout, 'Qout': Qout,
    'radius': radius, 'velo': velo, 'stress': stress, 'tau': tau, 'volume': volume,
    'dPin_dt': dPin_dt, 'dQin_dt': dQin_dt, 'dPout_dt': dPout_dt, 'dQout_dt': dQout_dt,
    'dradius_dt': dradius_dt, 'dvelo_dt': dvelo_dt, 'dstress_dt': dstress_dt,
    'dtau_dt': dtau_dt, 'dvolume_dt': dvolume_dt,
    'rho': rho, 'thick0': thick0, 'radius0': radius0,
    'C0': C0, 'C1': C1, 'C2': C2, 'C3': C3,
    'W1': W1, 'W2': W2,
    'eta': eta, 'sigma_max': sigma_max,
    'pi': sp.pi, 'M_PI': sp.pi, 'exp': sp.exp, 'pow': sp.Pow, 'fabs': sp.Abs,
    'sp': sp
}


def extract_cpp_expression(cpp_expr):
    """Convert C++ expression to Python/SymPy expression."""
    # Replace C++ specific syntax
    expr = cpp_expr.replace('pow(', 'Pow(')
    expr = expr.replace('M_PI', 'pi')
    expr = expr.replace('fabs(', 'Abs(')
    expr = expr.replace('std::abs(', 'Abs(')
    return expr


def parse_expression(expr_str):
    """Parse expression string to SymPy expression."""
    try:
        expr_clean = extract_cpp_expression(expr_str)
        expr = sp.sympify(expr_clean, locals=symbols_dict)
        return expr
    except Exception as e:
        print(f"Error parsing: {expr_str[:100]}...")
        print(f"Error: {e}")
        return None


def compare_equations(eq_growth, eq_expmat, eq_name):
    """Compare two equations and check if they're equivalent when Theta_c=1, Theta_r=1."""
    
    if eq_growth is None or eq_expmat is None:
        return None
    
    # Substitute Theta_c = 1 and Theta_r = 1
    eq_growth_sub = eq_growth.subs({Theta_c: 1, Theta_r: 1})
    
    # Simplify both sides
    eq_growth_simplified = simplify(eq_growth_sub)
    eq_expmat_simplified = simplify(eq_expmat)
    
    # Check if they're equivalent
    difference = simplify(eq_growth_simplified - eq_expmat_simplified)
    is_equivalent = (difference == 0)
    
    return {
        'name': eq_name,
        'equivalent': is_equivalent,
        'growth_original': eq_growth,
        'growth_after_sub': eq_growth_sub,
        'growth_simplified': eq_growth_simplified,
        'expmat_simplified': eq_expmat_simplified,
        'difference': difference
    }


def main():
    print("="*90)
    print("EQUATION EQUIVALENCE CHECK: ChamberSphere_growth vs ChamberSphere")
    print("="*90)
    print("Checking if ChamberSphere_growth equations reduce to ChamberSphere")
    print("when Theta_r = 1 and Theta_c = 1\n")
    
    # =========================================================================
    # EQUATION 0: Balance of Linear Momentum (Constant Coefficient)
    # =========================================================================
    print("\n" + "="*90)
    print("EQUATION 0: Balance of Linear Momentum - Constant Coefficient E[0,5]")
    print("="*90)
    
    coeff_growth_0 = parse_expression("pow(Theta_c, 2)*Theta_r*rho*thick0")
    coeff_expmat_0 = parse_expression("rho * thick0")
    
    result_0 = compare_equations(coeff_growth_0, coeff_expmat_0, "E[0,5]")
    print(f"\nGrowth:  {result_0['growth_original']}")
    print(f"Expmat:  {result_0['expmat_simplified']}")
    print(f"After substitution Theta_c=1, Theta_r=1:")
    print(f"  Growth becomes: {result_0['growth_simplified']}")
    print(f"  Difference: {result_0['difference']}")
    print(f"Result: {'✓ EQUIVALENT' if result_0['equivalent'] else '✗ NOT EQUIVALENT'}")
    
    # =========================================================================
    # EQUATION 1: Spherical Stress - Variable Coefficient C[1]
    # =========================================================================
    print("\n" + "="*90)
    print("EQUATION 1: Spherical Stress - Main Equation C[1]")
    print("="*90)
    
    # Growth version (extremely complex)
    eq_growth_1 = parse_expression(
        "2*(-4*C0*C1*(radius + radius0)*(2*pow(Theta_c, 2)*pow(radius + radius0, 6) + "
        "pow(Theta_r, 2)*pow(radius0, 6) - 3*pow(radius0, 2)*pow(radius + radius0, 4)*"
        "fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)))*(pow(Theta_r, 2)*"
        "pow(radius0, 6) - pow(Theta_r, 2)*pow(radius + radius0, 6) + (pow(Theta_c, 2)*pow(Theta_r, 2) - 1)*"
        "(0.66666666666666663*pow(Theta_c, 2)*pow(radius + radius0, 6) + "
        "0.33333333333333331*pow(Theta_r, 2)*pow(radius0, 6)))*exp(C1*pow(2*pow(Theta_c, 2)*"
        "pow(radius + radius0, 6) + pow(Theta_r, 2)*pow(radius0, 6) - 3*pow(radius0, 2)*"
        "pow(radius + radius0, 4)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)/("
        "pow(radius0, 4)*pow(radius + radius0, 8)*pow(fabs(pow(Theta_c, 1.3333333333333333)*"
        "pow(Theta_r, 0.66666666666666663)), 2))) + C2*C3*(1.3333333333333333*pow(Theta_c, 2) + "
        "0.66666666666666663*pow(Theta_r, 2))*pow(radius + radius0, 11)*(pow(Theta_c, 2)*"
        "pow(radius + radius0, 2) - pow(radius0, 2)*fabs(pow(Theta_c, 1.3333333333333333)*"
        "pow(Theta_r, 0.66666666666666663)))*exp(C3*pow(pow(Theta_c, 2)*pow(radius + radius0, 2) - "
        "pow(radius0, 2)*fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2)/("
        "pow(radius0, 4)*pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))) + "
        "pow(Theta_r, 2)*dradius_dt*eta*(pow(Theta_c, 2)*pow(radius + radius0, 12) + "
        "2*pow(Theta_r, 2)*pow(radius0, 12))*pow(fabs(pow(Theta_c, 1.3333333333333333)*"
        "pow(Theta_r, 0.66666666666666663)), 2))/(pow(Theta_r, 2)*pow(radius0, 2)*pow(radius + radius0, 11)*"
        "pow(fabs(pow(Theta_c, 1.3333333333333333)*pow(Theta_r, 0.66666666666666663)), 2))"
    )
    
    # Expmat version (simpler)
    eq_expmat_1 = parse_expression(
        "4 * (dradius_dt * eta * (-2 * pow(radius0, 12) + pow(radius + radius0, 12)) + "
        "pow(radius + radius0, 5) * (-pow(radius0, 6) + pow(radius + radius0, 6)) * "
        "(W1 * pow(radius0, 2) + W2 * pow(radius + radius0, 2))) / "
        "(pow(radius0, 2) * pow(radius + radius0, 11))"
    )
    
    # Simplify at Theta_c=1, Theta_r=1
    print("\nNote: This is the most complex equation. Simplifying...")
    if eq_growth_1 is not None:
        eq_growth_1_sub = eq_growth_1.subs({Theta_c: 1, Theta_r: 1})
        eq_growth_1_simp = simplify(eq_growth_1_sub)
        print(f"\nGrowth equation (after Theta_c=1, Theta_r=1 substitution, simplified):")
        print(f"  {eq_growth_1_simp}")
    
    if eq_expmat_1 is not None:
        eq_expmat_1_simp = simplify(eq_expmat_1)
        print(f"\nExpmat equation (simplified):")
        print(f"  {eq_expmat_1_simp}")
    
    if eq_growth_1 is not None and eq_expmat_1 is not None:
        result_1 = compare_equations(eq_growth_1, eq_expmat_1, "C[1]")
        print(f"\nEquivalence check:")
        print(f"  Difference after simplification: {result_1['difference']}")
        print(f"Result: {'✓ EQUIVALENT' if result_1['equivalent'] else '✗ NOT EQUIVALENT'}")
    else:
        result_1 = None
        print("Could not parse equations for full comparison.")
    
    # =========================================================================
    # EQUATION 2: Volume Change - Main Equation C[2]
    # =========================================================================
    print("\n" + "="*90)
    print("EQUATION 2: Volume Change - Main Equation C[2]")
    print("="*90)
    
    eq_growth_2 = parse_expression("4*M_PI*velo*pow(radius + radius0, 2)")
    eq_expmat_2 = parse_expression("4 * M_PI * velo * pow(radius + radius0, 2)")
    
    result_2 = compare_equations(eq_growth_2, eq_expmat_2, "C[2]")
    print(f"\nGrowth:  {result_2['growth_original']}")
    print(f"Expmat:  {result_2['expmat_simplified']}")
    print(f"Result: {'✓ EQUIVALENT' if result_2['equivalent'] else '✗ NOT EQUIVALENT'}")
    
    # =========================================================================
    # EQUATION 3: Active Stress - Main Equation C[3]
    # =========================================================================
    print("\n" + "="*90)
    print("EQUATION 3: Active Stress - Main Equation C[3]")
    print("="*90)
    
    eq_growth_3 = parse_expression("-act_plus*sigma_max")
    eq_expmat_3 = parse_expression("-act_plus * sigma_max")
    
    result_3 = compare_equations(eq_growth_3, eq_expmat_3, "C[3]")
    print(f"\nGrowth:  {result_3['growth_original']}")
    print(f"Expmat:  {result_3['expmat_simplified']}")
    print(f"Result: {'✓ EQUIVALENT' if result_3['equivalent'] else '✗ NOT EQUIVALENT'}")
    
    # =========================================================================
    # EQUATION 0: Balance of Momentum - Main Equation C[0]
    # =========================================================================
    print("\n" + "="*90)
    print("EQUATION 0: Balance of Momentum - Main Equation C[0]")
    print("="*90)
    
    # Growth version
    eq_growth_0_main = parse_expression(
        "Theta_c*(radius + radius0)*(-Pout*Theta_c*(radius + radius0) + stress*thick0)/pow(radius0, 2)"
    )
    
    # Expmat version
    eq_expmat_0_main = parse_expression(
        "(radius + radius0) * (-Pout * (radius + radius0) + stress * thick0) / pow(radius0, 2)"
    )
    
    result_0_main = compare_equations(eq_growth_0_main, eq_expmat_0_main, "C[0]")
    print(f"\nGrowth:  {result_0_main['growth_original']}")
    print(f"Expmat:  {result_0_main['expmat_simplified']}")
    print(f"After substitution Theta_c=1, Theta_r=1:")
    print(f"  Growth becomes: {result_0_main['growth_simplified']}")
    print(f"  Difference: {result_0_main['difference']}")
    print(f"Result: {'✓ EQUIVALENT' if result_0_main['equivalent'] else '✗ NOT EQUIVALENT'}")
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("\n" + "="*90)
    print("SUMMARY")
    print("="*90)
    
    results = [result_0, result_2, result_3, result_0_main]
    if result_1 is not None:
        results.insert(1, result_1)
    
    equiv_count = sum(1 for r in results if r and r['equivalent'])
    total_count = len([r for r in results if r is not None])
    
    print(f"\nEquations checked: {total_count}")
    print(f"Equivalent equations: {equiv_count}")
    print(f"Non-equivalent equations: {total_count - equiv_count}\n")
    
    for r in results:
        if r:
            status = "✓ EQUIVALENT" if r['equivalent'] else "✗ NOT EQUIVALENT"
            print(f"  {r['name']}: {status}")
    
    if equiv_count == total_count and total_count > 0:
        print("\n" + "="*90)
        print("✓ SUCCESS!")
        print("="*90)
        print("All checked equations are mathematically equivalent!")
        print("ChamberSphere_growth with Theta_r=1 and Theta_c=1 reduces to")
        print("ChamberSphere (the expmat version without growth).")
    elif total_count == 0:
        print("\n✗ ERROR: No equations were successfully parsed")
    else:
        print("\n" + "="*90)
        print("⚠ PARTIAL EQUIVALENCE")
        print("="*90)
        print("Some equations are NOT equivalent. Review the differences above.")


if __name__ == "__main__":
    main()