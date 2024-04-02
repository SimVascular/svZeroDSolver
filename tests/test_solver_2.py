from .conftest import run_with_reference
from .utils import execute_pysvzerod

def test_steadyFlow_R_R():
    
    ref, config = execute_pysvzerod("tests/cases/steadyFlow_R_R.json", "solver")

    run_with_reference(ref, config)


if __name__ == "__main__":
    test_steadyFlow_R_R()
