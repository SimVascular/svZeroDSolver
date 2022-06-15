# coding=utf-8

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

import numpy as np
from scipy.sparse import linalg


class GeneralizedAlpha:
    """Generalized alpha time integrator.

    Solves system E*ydot + F*y + C = 0 with generalized alpha and Newton-Raphson for non-linear residual
    """

    def __init__(self, rho, n, time_step_size):

        # Setup constants for generalized alpha time integration
        self.alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho)
        self.alpha_f = 1.0 / (1.0 + rho)
        self.gamma = 0.5 + self.alpha_m - self.alpha_f

        self.fac = self.alpha_m / (self.alpha_f * self.gamma)

        # problem dimension
        self.n = n
        self.time_step_size = time_step_size
        self.time_step_size_inv = 1.0 / time_step_size

        # jacobian matrix
        if self.n > 800:
            self.solve = linalg.spsolve
        else:
            self.solve = np.linalg.solve

        # stores matrices E, F, vector C, and tangent matrices dE, dF, dC
        self.mat = {}
        self.mats = ["E", "F", "dE", "dF", "dC"]
        self.vecs = ["C"]
        for m in self.mats:
            self.mat[m] = np.zeros((self.n, self.n))
        for v in self.vecs:
            self.mat[v] = np.zeros(self.n)

    def assemble(self, block_list, y, ydot):
        """
        Assemble block matrices into global matrices
        """
        for bl in block_list:
            while bl.vecs_to_assemble:
                n = bl.vecs_to_assemble.pop()
                self.mat[n][bl.global_row_id] = bl.vec[n]
            while bl.mats_to_assemble:
                n = bl.mats_to_assemble.pop()
                self.mat[n][bl.flat_row_ids, bl.flat_col_ids] = bl.mat[n].ravel()

        lhs = self.mat["F"] + (
            self.mat["dE"]
            + self.mat["dF"]
            + self.mat["dC"]
            + self.mat["E"] * self.fac * self.time_step_size_inv
        )
        rhs = -self.mat["E"].dot(ydot) - self.mat["F"].dot(y) - self.mat["C"]

        return lhs, rhs

    def step(self, y, ydot, time, block_list, max_iter=30):
        """
        Perform one time step
        """

        # Make copies to prevent overwriting
        y = y.copy()
        ydot = ydot.copy()

        # initial guess for time step
        curr_y = y + 0.5 * self.time_step_size * ydot
        curr_ydot = ydot * ((self.gamma - 0.5) / self.gamma)

        # Substep level quantities
        yaf = y + self.alpha_f * (curr_y - y)
        ydotam = ydot + self.alpha_m * (curr_ydot - ydot)

        # initialize solution
        time = time + self.alpha_f * self.time_step_size

        # initialize blocks
        for b in block_list:
            b.update_time(time)

        fac_ydotam = self.fac * self.time_step_size_inv
        for iter in range(max_iter):
            # update solution-dependent blocks
            for b in block_list:
                b.update_solution(yaf)

            # update residual and jacobian
            lhs, res = self.assemble(block_list, yaf, ydotam)

            # solve for Newton increment
            dy = self.solve(lhs, res)

            # update solution
            yaf += dy
            ydotam += dy * fac_ydotam

            if np.abs(res).max() <= 5e-4:
                break

        if iter == max_iter - 1:
            print(
                f"Max NR iterations reached at time: {time:.3f}s with, max error: {max(abs(res))}"
            )

        # update time step
        curr_y = y + (yaf - y) / self.alpha_f
        curr_ydot = ydot + (ydotam - ydot) / self.alpha_m

        return curr_y, curr_ydot
