# Copyright (c) Stanford University, The Regents of the University of
#               California, and others.
#
# All Rights Reserved.
#
# See Copyright-SimVascular.txt for additional details.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""This module holds the DOFHandler class."""


class DOFHandler:
    """Degree-of-freedom handler.

    This class handles degrees-of-freedom for model variables and equations.
    It assigns each element with row and column indices which it can use to
    assemble it's local contributions into the local contributions into the
    global system.

    Attributes:
        variables: List of variable names corresponding to the global IDs.
            Variables without a name have an entry None.
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
