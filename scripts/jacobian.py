from sympy import symbols, Matrix, diff, simplify, pi
from sympy.printing import ccode
import re
import pdb

# --- Global Variable Definitions ---
ny = 9
y = Matrix(symbols(f'y0:{ny}'))
dy = Matrix(symbols(f'dy0:{ny}'))

# Parameter symbols
rho, d, Ro, W1, W2, eta, a, a_plus, sigma_o = symbols('rho d Ro W1 W2 eta a a_plus sigma_o')
time_dependent_symbols = {a, a_plus}

# Helper functions
def C(Ro, idx):
    return (1 + (y[idx] / Ro)) ** 2

def dC(Ro, idx, d_idx):
    return diff(C(Ro, idx), y[idx]) * dy[d_idx]

# --- User-Defined Residuals ---
residuals = Matrix([
    rho * d * dy[5] + (d / Ro) * (1 + (y[4] / Ro)) * y[6] - y[2] * C(Ro, 4),
    -y[6] + y[7] + 4 * (1 - C(Ro, 4) ** -3) * (W1 + C(Ro, 4) * W2) + 2 * eta * dC(Ro, 4, 4) * (1 - 2 * C(Ro, 4) ** -6),
    dy[7] + a * y[7] - sigma_o * a_plus,
    y[1] - y[3] - dy[8],
    4 * pi * Ro ** 2 * C(Ro, 4) * y[5] - dy[8],
    dy[4] - y[5],
    y[0] - y[2]
])

#
nr = residuals.shape[0]
assert ny - 2 == nr, f"Number of residuals ({nr}) should be number of unknowns ({ny}) minus 2."

def extract_linear_terms(residuals, y, dy):
    E = Matrix.zeros(nr, ny)
    F = Matrix.zeros(nr, ny)
    def is_constant_coeff(coeff):
        return (coeff and coeff.free_symbols.isdisjoint(y.free_symbols | dy.free_symbols))
    for i in range(nr):
        for j in range(ny):
            for mat, dd in zip([E, F], [dy, y]):
                coeff = residuals[i].coeff(dd[j])
                if is_constant_coeff(coeff):
                    mat[i, j] = coeff
    return E, F

def compute_nonlinear_part(residuals, E, F, y, dy):
    c = simplify(residuals - E * dy - F * y)
    dc_dy = simplify(c.jacobian(y))
    dc_ddy = simplify(c.jacobian(dy))
    return c, dc_dy, dc_ddy

def partition_terms(E, F, c, dc_dy, dc_ddy, time_dependent_symbols):
    def depends_on(expr, symbols_set):
        return any(sym in expr.free_symbols for sym in symbols_set)
    parts = {"constant": [], "time_dependent": [], "solution_dependent": []}
    for label, mat in [("E", E), ("F", F)]:
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if mat[i, j] != 0:
                    target = "time_dependent" if depends_on(mat[i, j], time_dependent_symbols) else "constant"
                    parts[target].append((label, i, j, mat[i, j]))
    for i in range(nr):
        if c[i] != 0:
            parts["solution_dependent"].append(("c", i, c[i]))
    return parts

def replace_symbolic_indices(expr, base):
    pattern = re.compile(rf'\b{base}(\d+)\b')
    return pattern.sub(rf'{base}[global_var_ids[\1]]', expr)

def format_cpp_expr(expr):
    cpp_expr = ccode(expr).replace('3.141592653589793', 'M_PI')
    cpp_expr = replace_symbolic_indices(cpp_expr, 'y')
    cpp_expr = replace_symbolic_indices(cpp_expr, 'dy')
    return cpp_expr

def print_cpp(system_part, i, j=None, expr=None):
    if expr != 0:
        if j is not None:
            print(f"  system.{system_part}[{i}][{j}] = {format_cpp_expr(expr)};")
        else:
            print(f"  system.{system_part}[{i}] = {format_cpp_expr(expr)};")

def generate_cpp_function(parts):
    for part in parts:
        label, i, j_or_expr, expr_or_none = part if len(part) == 4 else (part[0], part[1], None, part[2])
        if label in ["E", "F"]:
            print_cpp(label, i, j_or_expr, expr_or_none)
        elif label == "c":
            print_cpp('c', i, expr=j_or_expr)

def generate_cpp_solution_function(c, dc_dy, dc_ddy):
    for i in range(nr):
        print_cpp('c', i, expr=c[i])
        for j in range(ny):
            for mat, name in zip([dc_dy, dc_ddy], ['dc_dy', 'dc_ddy']):
                print_cpp(name, i, j, mat[i, j])

def main():
    E, F = extract_linear_terms(residuals, y, dy)
    c, dc_dy, dc_ddy = compute_nonlinear_part(residuals, E, F, y, dy)
    parts = partition_terms(E, F, c, dc_dy, dc_ddy, time_dependent_symbols)
    print('update_constant')
    generate_cpp_function(parts['constant'])
    print('update_time')
    generate_cpp_function(parts['time_dependent'])
    print('update_solution')
    generate_cpp_solution_function(c, dc_dy, dc_ddy)

main()