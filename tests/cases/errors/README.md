# Input files that fail to load

Each file here is a small, otherwise valid model with exactly one mistake in it.
They are examples of the errors svZeroDSolver reports for a model that is not
connected the way its blocks require, and reproducers for the checks that report
them.

Run one with:

```bash
svzerodsolver tests/cases/errors/vessel_not_connected.json out.csv
```

Every file carries its mistake and the error it produces under `description`.
The solver reports the error and exits with code 1; none of them run.

| File | Mistake |
| --- | --- |
| `vessel_not_connected.json` | A vessel in no junction and with no boundary conditions |
| `vessel_connected_on_one_side.json` | A vessel with a boundary condition on one side only |
| `boundary_condition_not_connected.json` | A boundary condition no vessel names |
| `junction_unknown_vessel_id.json` | A junction listing a `vessel_id` no vessel defines |
| `junction_missing_connection_keys.json` | A junction giving `inlet_vessels` but no `outlet_vessels` |

The first three are caught by the connection count check in `Block::setup_dofs_`:
a block contributes one equation per variable it owns, which is one per
connection plus one per internal variable, so a block whose connections do not
add up cannot assemble its equations. Before that check it indexed past the end
of `global_var_ids` and the solver segfaulted.
