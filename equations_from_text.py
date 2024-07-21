import sympy as sp


def read_equations_from_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    equations = []
    for line in lines:
        # Remove any whitespace or newline characters
        equation = line.strip()
        if equation:
            # Check if the equation contains an '=' for assignment
            if '=' in equation:
                lhs, rhs = equation.split('=')
                lhs = lhs.strip()
                rhs = rhs.strip()
                # Create a sympy symbol for the left-hand side
                lhs_symbol = sp.symbols(lhs)
                # Convert the right-hand side to a sympy expression
                rhs_expr = sp.sympify(rhs)
                # Form the equation
                sympy_eq = sp.Eq(lhs_symbol, rhs_expr)
            else:
                # Convert the string to a sympy equation
                sympy_eq = sp.sympify(equation)
            equations.append(sympy_eq)
    
    return equations


if __name__ == "__main__":

    # Example usage
    file_path = './docs/equations.txt'
    equations = read_equations_from_file(file_path)

    for eq in equations:
        print(eq.lhs, '=', eq.rhs)
        #print(eq.free_symbols)
