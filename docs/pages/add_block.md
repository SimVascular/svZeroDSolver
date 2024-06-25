@page add_block Adding New Blocks

[TOC]

Below are details on the steps required to implement a new block in svZeroDSolver.

*Note: The best way to implement a new block is to look at examples of existing block classes. See the `ValveTanh` class for an example.*

# 1. Add the new block to the relevant lists/dictionaries.

* `BlockType` in src/model/BlockType.h
* `block_factory_map` in src/model/Model.cpp
  * *Note: In `block_factory_map`, the dictionary key should match the string specifying the type of block in the `.json` configuration/input file, and the dictionary value should match the class constructor name for the block.*
* If the new block requires special handling that is different from the current blocks (most new blocks do not), add a new category to `BlockClass` in src/model/BlockType.h

<p> <br> </p>

# 2. Create a class for the new block 

* The new class will be inherited from `Block`. Define a constructor of the form:
```
    GenericBlock(int id, Model *model) 
    : Block(id, model, BlockType::block_type, BlockClass::block_class,
    {{Param_1, InputParameter()}, 
    {Param_2, InputParameter()}, 
    ..., 
    {Param_N, InputParameter()}}) {}
```
  * `GenericBlock` is the name of the new class
  * `block_type` and `block_class` are the same as what was added in Step 1 above.
  * The names of the input parameters of the block are `Param_1`, ... , `Param_N`. 
  * The properties of each parameter are defined by `InputParameter`, which specifies whether it is optional, an array, a scalar, a function, and its default value.
  * The names `Param_1`, ... , `Param_N` should be the same as the parameter names within the block definition in the `.json` configuration/input file.

<p> <br> </p>

* The class should have a `setup_dofs(DOFHandler &dofhandler)` function.
  * This function typically only includes a call to the following function:
  ```
  Block::setup_dofs_(DOFHandler &dofhandler, int num_equations, const std::list<std::string> &internal_var_names)
  ```
  * In the above function, `num_equations` is the number of governing equations for the new block. 
  * `internal_var_names` is a list of strings that specify names for variables that are internal to the block, i.e. all variables for the block apart from the flow and pressure at the block's inlets and outlets.

<p> <br> </p>

* The class should have a `TripletsContributions num_triplets{*, *, *}` object. 
  * This specifies how many elements the governing equations of the block contribute to the global `F`, `E` and `dC_dy` matrices respectively. Details are in Step 3 below. 

<p> <br> </p>

* The class should have an `update_constant` function and may also contain `update_time` and `update_solution` functions. These functions implement the governing equations for the block. Details are in Steps 3-4 below.

<p> <br> </p>

* *(Optional)* The class can have an  `enum ParamId` object that relates the parameter indices to their names. 
  * This makes it easier to reference the parameters while implementing the governing equations of the block (discussed below). 
  * The order of parameters in the `ParamId` object should match the order in the constructor. 

<p> <br> </p>

# 3. Set up the governing equations for the block.

* The local state vector for each block is always arranged as `[P_in, Q_in, P_out, Q_out, InternalVariable_1, ..., InternalVariable_N]`. 
  
  * Here, `InternalVariable*` refers to any variable in the governing equations that are not the inlet and outlet flow and pressure. These are the same as those discussed above in the function `setup_dofs`.
  
  * The length of the state vector is typically four (inlet and outlet pressure and flow) plus the number of internal variables.

<p> <br> </p>

* The equations should be written in the form `E(t)*ydot + F(t)*y + C(y,t) = 0`. 
  
  * `y` is the local state vector mentioned above 
  
  * `ydot` is the time-derivative of the local state vector 
  
  * `E` and `F` are matrices of size `number_of_equations*size_of_state_vector`. 
  
  * `c` is a vector of length `number_of_equations`. 
  
  * `E` and `F` contain terms of the governing equation that multiply the respective components of `ydot` and `y` respectively.
  
  * `C` is a vector containing all non-linear and constant terms in the equation. 
  
  * If the equation contains non-linear terms, the developer should also derive the derivative of `C` with respect to `y` and `ydot`. These will be stored in the block's `dC_dy` and `dC_dydot` matrices, both of which are size `number_of_equations*size_of_state_vector`.

<p> <br> </p>

* For example, let's assume a block has the following non-linear governing equations:
```
a*dQ_in/dt + b*P_in + c*(dP_in/dt)*Q_in + d = 0
```
```
e*dP_out/dt + f*Q_out*Q_out + g*P_out + h*I_1 = 0
```
  
  * For this block, the `P_in` and `Q_in` are the pressure and flow at the inlet respectively, `P_out` and `Q_out` are the pressure and flow at the outlet, and `I_1` is an internal variable.
  
  * The state vector is `[P_in, Q_in, P_out, Q_out, I_1]`.
  
  * The contributions to the local `F` matrix are `F[0,0] = b`, `F[1,2] = g` and `F[1,4] = h`.
  
  * The contributions to the local `E` matrix are `E[0,1] = a` and `E[1,2] = e`.
  
  * The contributions to the local `C` vector are `C[0] = c*(dP_in/dt)*Q_in + d` and `C[1] = f*Q_out*Q_out`.
  
  * The contributions to the local `dC_dy` matrix are `dC_dy[0,1] = c*(dP_in/dt)` and `dC_dy[1,3] = 2*f*Q_out`.
  
  * The contributions to the local `dC_dydot` matrix are `dC_dydot[0,0] = c*Q_in`.
  
  * In this case, the block has 3 contributions to `F`, 2 contributions to `E`, and 2 constributions to `dC_dy`. So the class will have a member `TripletsContributions num_triplets{3, 2, 2}`.

# 4. Implement the matrix equations for the block.

* Implement the `update_constant`, `update_time` and `update_solution` functions.

  * All matrix elements that are constant are specified in `update_constant`.

  * Matrix elements that depend only on time (not the solution of the problem itself) are specified in `update_time`. 

  * Matrix elements that change with the solution (i.e. depend on the state variables themselves) are specified in `update_solution`. 
  
  * Not all blocks will require the `update_time` and `update_solution` functions.

<p> <br> </p>

* The elements of the `E`, `F`, `dC_dy` and `dC_dydot` matrices are populated using the syntax 
```
system.F.coeffRef(global_eqn_ids[current_block_equation_id], global_var_ids[current_block_variable_ids]) = a
```
  
  * Here, `current_block_equation_id` goes from 0 to `number_of_equations-1` (for the current block) and `current_block_variable_ids` goes from 0 to `size_of_state_vector-1` for the current block. 

<p> <br> </p>

* If the block contains non-linear equations, these terms must be specified in `update_solution` as `system.C(global_eqn_ids[current_block_equation_id]) = non_linear_term`.

<p> <br> </p>

* For non-linear equations, the derivative of the term specified above with respect to each state variable must also be provided. This goes into a `dC_dy` matrix using the following syntax 
```
system.dC_dy.coeffRef(global_eqn_ids[current_block_equation_id], global_var_ids[current_block_variable_id]) = a
```
  
  * Here, `a` is the derivative of the non-linear term in the equation with ID `current_block_equation_id` with respect to the local state variable with ID `current_block_variable_id`. 
  
  * For example, if the non-linear term is in the first equation, `current_block_equation_id = 0`. 
  
  * For the derivative of this term with respect to `P_in`, `current_block_variable_id = 0` and for the derivative of this term with respect to `P_out`, `current_block_variable_id = 2`.

<p> <br> </p>

# 4. Add the new block to the build system 

* Add `GenericBlock.h` and `GenericBlock.cpp` to `src/model/CMakeLists.txt`

