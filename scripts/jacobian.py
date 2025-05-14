from sympy import symbols, Matrix, simplify, pi
from sympy.printing import ccode
import re
import pdb
import yaml

def load_model(filepath):
    with open(filepath, 'r') as file:
        data = yaml.safe_load(file)

    variables = Matrix(symbols(' '.join(data['variables'])))
    derivatives = Matrix(symbols(' '.join(data['derivatives'])))
    constants = symbols(' '.join(data['constants']))

    context = {str(s): s for s in list(variables) + list(derivatives) + list(constants)}
    context['pi'] = pi

    exec(data['helper_functions'], context, context)
    residual_exprs = [eval(res, context, context) for res in data['residuals']]

    residuals = Matrix(residual_exprs)
    time_dependent = {context[s] for s in data.get('time_dependent', [])}

    assert time_dependent.issubset({s for s in constants}), "Time dependent variables must be constants"
    assert len(variables) == len(derivatives), f"Number of variables must be equal to number of derivatives"
    assert len(variables) - 2 == len(residuals), f"Number of residuals must be number of unknowns minus 2"

    return variables, derivatives, constants, residuals, time_dependent

def extract_linear(residuals, y, dy):
    def is_constant_coeff(coeff):
        return (coeff and coeff.free_symbols.isdisjoint(y.free_symbols | dy.free_symbols))
    nr = residuals.shape[0]
    ny = y.shape[0]
    E = Matrix.zeros(nr, ny)
    F = Matrix.zeros(nr, ny)
    for i in range(nr):
        for j in range(ny):
            for mat, dd in zip([E, F], [dy, y]):
                coeff = residuals[i].coeff(dd[j])
                if is_constant_coeff(coeff):
                    mat[i, j] = coeff
    return E, F

def extract_nonlinear(residuals, E, F, y, dy):
    c = simplify(residuals - E * dy - F * y)
    dc_dy = simplify(c.jacobian(y))
    dc_ddy = simplify(c.jacobian(dy))
    return c, dc_dy, dc_ddy

def depends_on(expr, symbols_set):
    return any(sym in expr.free_symbols for sym in symbols_set)

def partition_terms(E, F, c, dc_dy, dc_ddy, time_dependent_symbols):
    parts = {"constant": [], "time": [], "solution": []}
    for i in range(E.shape[0]):
        parts["solution"].append(("c", i, c[i]))
        for j in range(E.shape[1]):
            for label, mat in [("E", E), ("F", F), ("dc_dy", dc_dy), ("dc_ddy", dc_ddy)]:
                if label.startswith('d'):
                    target = "solution"
                elif depends_on(mat[i, j], time_dependent_symbols):
                    target = "time"
                else:
                    target = "constant"
                parts[target].append((label, i, j, mat[i, j]))
    return parts

def replace_symbolic_indices(expr, base):
    pattern = re.compile(rf'\b{base}(\d+)\b')
    return pattern.sub(rf'{base}[global_var_ids[\1]]', expr)

def format_cpp_expr(expr):
    cpp_expr = ccode(expr).replace('3.141592653589793', 'M_PI')
    cpp_expr = replace_symbolic_indices(cpp_expr, 'y')
    cpp_expr = replace_symbolic_indices(cpp_expr, 'dy')
    return cpp_expr

def print_system(parts):
    for part in parts:
        if part[-1] != 0:
            str = f"  system.{part[0]}[{part[1]}]"
            if len(part) == 4:
                str += f"[{part[2]}]"
            print(str + f" = {format_cpp_expr(part[-1])};")

def print_constants(parts):

    pdb.set_trace()

def main(yaml_path):
    y, dy, constants, residuals, time_dependent = load_model(yaml_path)
    E, F = extract_linear(residuals, y, dy)
    c, dc_dy, dc_ddy = extract_nonlinear(residuals, E, F, y, dy)
    parts = partition_terms(E, F, c, dc_dy, dc_ddy, time_dependent)
    for section in ['constant', 'time', 'solution']:
        print(f'update_{section}')
        print_system(parts[section])
        print('\n')

main("ChamberSphere.yaml")