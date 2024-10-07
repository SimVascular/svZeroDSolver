import pytest
import filecmp
import os
import sys

# Append the 'applications/dirgraph' directory to the system path for imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../applications/svZeroDVisualization')))

# Import the function from applications/dirgraph
from dirgraph_utils import set_up_0d_network


cases_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'cases'))
results_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'cases/dirgraph-results'))


# List of files to exclude
excluded_files = [
    'closedLoopHeart_singleVessel.json',
    'closedLoopHeart_withCoronaries.json',
    'coupledBlock_closedLoopHeart_singleVessel.json',
    'coupledBlock_closedLoopHeart_withCoronaries.json'
]

# Generate the list of JSON files to test
json_files = [
    filename for filename in os.listdir(cases_dir)
    if filename.endswith('.json') and filename not in excluded_files
]

@pytest.fixture
def setup_files(tmp_path, filename):
    # Setup input file path
    input_file = os.path.join(cases_dir, filename)
    expected_output_file = results_dir + '/' + os.path.splitext(filename)[0] + "_directed_graph.dot"
    return input_file, expected_output_file, tmp_path


def run_program(zero_d_solver_input_file_path, tmp_path):
    set_up_0d_network(
        zero_d_solver_input_file_path,
        name_type='id',
        draw_directed_graph=False,
        output_dir = tmp_path
    )


@pytest.mark.parametrize("filename", json_files)
def test_directed_graph_generation(setup_files):
    input_file_path, expected_dot_file_path, tmp_path = setup_files

    # Run your program with the input file
    run_program(input_file_path, tmp_path)

    generated_dot_file_path = tmp_path / (os.path.splitext(os.path.basename(input_file_path))[0] + "_directed_graph.dot")

    filecmp.clear_cache()
    assert filecmp.cmp(generated_dot_file_path, expected_dot_file_path), \
        f"The generated dot file '{generated_dot_file_path}' does not match the expected dot file '{expected_dot_file_path}'."

if __name__ == "__main__":
    pytest.main()
