class DOFHandler:
    """Degree-of-freedom handler.

    Class for handling global degrees of freedom of a system. Assigns global
    IDs to variables and equations.

    Attributes:
        variables: List of variable names corresponding to the global IDs.
            Variables without a name have an entry None.
        n: Size of the system.
    """

    def __init__(self):
        """Create a new DOFHandler instance."""
        self._var_counter = -1
        self._eqn_counter = -1
        self.variables: list[str] = []

    @property
    def n(self) -> int:
        """Size of the system."""
        return self._eqn_counter + 1

    def register_variable(self, name: str = None) -> int:
        """Register a new variable and get global index.

        Args:
            name: Optional name of the variables.

        Returns:
            global_id: Global id of the variable.
        """
        self._var_counter += 1
        self.variables.append(name)
        return self._var_counter

    def register_equation(self) -> int:
        """Register a new equation and get global index.

        Returns:
            global_id: Global id of the equation.
        """
        self._eqn_counter += 1
        return self._eqn_counter
