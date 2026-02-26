from sympy import symbols, Matrix, simplify, pi, Abs
from sympy.printing import ccode
import re
import pdb
import yaml
import argparse

def load_model(filepath):
    with open(filepath, 'r') as file:
        data = yaml.safe_load(file)

    variables = Matrix(symbols(' '.join(data['variables']), real=True))
    derivatives = Matrix(symbols(' '.join(data['derivatives']), real=True))
    constants = symbols(' '.join(data['constants']), real=True)

    context = {str(s): s for s in list(variables) + list(derivatives) + list(constants)}
    context['pi'] = pi
    context['abs'] = Abs

    if 'helper_functions' in data:
        exec(data['helper_functions'], context, context)
    residual_exprs = [eval(res, context, context) for res in data['residuals']]

    residuals = Matrix(residual_exprs)
    time_dependent = {context[s] for s in data.get('time_dependent', [])}

    assert len(variables) == len(derivatives), f"Number of variables must be equal to number of derivatives"

    return variables, derivatives, constants, residuals, time_dependent

def extract_linear(residuals, y, dy):
    zero_subs = [(v, 0) for v in list(y) + list(dy)]
    nr = residuals.shape[0]
    ny = y.shape[0]
    E = Matrix.zeros(nr, ny)
    F = Matrix.zeros(nr, ny)
    for i in range(nr):
        for j in range(ny):
            for mat, dd in zip([E, F], [dy, y]):
                coeff = residuals[i].coeff(dd[j])
                const_coeff = simplify(coeff.subs(zero_subs))
                if const_coeff:
                    mat[i, j] = const_coeff
    return E, F

def extract_nonlinear(residuals, E, F, y, dy):
    C = simplify(residuals - E * dy - F * y)
    dC_dy = simplify(C.jacobian(y))
    dC_dydot = simplify(C.jacobian(dy))
    return C, dC_dy, dC_dydot

def depends_on(expr, symbols_set):
    return any(sym in expr.free_symbols for sym in symbols_set)

def partition_terms(E, F, C, dC_dy, dC_dydot, variables, time_dependent_symbols):
    parts = {"constant": [], "time": [], "solution": []}
    variable_symbols = variables.free_symbols
    for label, mat in [("E", E), ("F", F), ("C", C), ("dC_dy", dC_dy), ("dC_dydot", dC_dydot)]:
        for i in range(mat.rows):
            for j in range(mat.cols):
                if mat.cols == 1:
                    expr = mat[i]
                else:
                    expr = mat[i, j]
                if expr != 0:
                    if depends_on(expr, variable_symbols):
                        target = "solution"
                    elif depends_on(expr, time_dependent_symbols):
                        target = "time"
                    else:
                        target = "constant"
                    if mat.cols == 1:
                        parts[target].append((label, i, expr))
                    else:
                        parts[target].append((label, i, j, expr))
    return parts

def get_expressions(parts, type):
    expressions = [part[-1] for part in parts]
    deps = set().union(*(expr.free_symbols for expr in expressions))
    return deps.intersection(type)

def replace_symbolic_indices(expr, base):
    pattern = re.compile(rf'\b{base}(\d+)\b')
    return pattern.sub(rf'{base}[global_var_ids[\1]]', expr)

def format_cpp_expr(expr):
    try:
        cpp_expr = ccode(expr).replace('3.141592653589793', 'M_PI')
    except:
        pdb.set_trace()
    cpp_expr = replace_symbolic_indices(cpp_expr, 'y')
    cpp_expr = replace_symbolic_indices(cpp_expr, 'dy')
    return cpp_expr

def print_index(i, j=None):
    print(f".coeffRef(global_eqn_ids[{i}]", end='')
    if j is not None:
        print(f", global_var_ids[{j}]", end='')
    print(f")", end='')

def print_system(parts):
    for part in parts:
        if part[-1] != 0:
            print(f"  system.{part[0]}", end='')
            if len(part) == 3:
                print_index(part[1])
            elif len(part) == 4:
                print_index(part[1], part[2])
            print(f" = {format_cpp_expr(part[-1])};")

def print_constants(parts, constants, time_dependent):
    for out in get_expressions(parts, constants):
        if out in time_dependent:
            print(f"  // compute time dependent constant {out}")
        else:
            print(f"  const double {out} = parameters[global_param_ids[ParamId::{out}]];")

def print_variables(parts, y, dy):
    for out in get_expressions(parts, y.free_symbols | dy.free_symbols):
        for name, vec in zip(['y', 'dy'], [y, dy]):
            if out in vec:
                index = next(i for i, sym in enumerate(vec) if sym == out)
                print(f"  const double {out} = {name}[global_var_ids[{index}]];")

def main(yaml_path):
    # read model from yaml file
    y, dy, constants, residuals, time_dependent = load_model(yaml_path)

    # extract linear and nonlinear terms
    E, F = extract_linear(residuals, y, dy)
    C, dC_dy, dC_dydot = extract_nonlinear(residuals, E, F, y, dy)

    # split into constant, time dependent and solution parts
    parts = partition_terms(E, F, C, dC_dy, dC_dydot, Matrix(list(y) + list(dy)), time_dependent)

    # print C++ code
    for section in ['constant', 'time', 'solution']:
        print(f'update_{section}')
        print_constants(parts[section], constants, time_dependent)
        print_variables(parts[section], y, dy)
        print_system(parts[section])
        print()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process YAML model file to generate C++ code for svZeroDSolver')
    parser.add_argument('yaml_file', help='Path to YAML model file')
    args = parser.parse_args()
    main(args.yaml_file)