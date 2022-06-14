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
import scipy
import scipy.sparse.linalg
import copy


class GenAlpha:
    """
    Solves system E*ydot + F*y + C = 0 with generalized alpha and Newton-Raphson for non-linear residual
    """

    def __init__(self, rho, y):
        # Constants for generalized alpha
        self.alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho)
        self.alpha_f = 1.0 / (1.0 + rho)
        self.gamma = 0.5 + self.alpha_m - self.alpha_f

        self.fac = self.alpha_m / (self.alpha_f * self.gamma)

        # problem dimension
        self.n = y.shape[0]

        # stores matrices E, F, vector C, and tangent matrices dE, dF, dC
        self.mat = {}

        # jacobian matrix
        self.M = np.zeros((self.n, self.n))
        self.sparse = False
        if self.n > 800:
            self.solver = scipy.sparse.linalg.spsolve
        else:
            self.solver = np.linalg.solve

        # residual vector
        self.res = np.zeros(self.n)

        self.mats = ["E", "F", "dE", "dF", "dC"]
        self.vecs = ["C"]
        for m in self.mats:
            self.mat[m] = np.zeros((self.n, self.n))
        for v in self.vecs:
            self.mat[v] = np.zeros(self.n)

    def assemble_structures(self, block_list):
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

    def form_matrix_NR(self, invdt):
        """
        Create Jacobian matrix
        """
        self.M = self.mat["F"] + (
            self.mat["dE"]
            + self.mat["dF"]
            + self.mat["dC"]
            + self.mat["E"] * self.fac * invdt
        )

    def form_rhs_NR(self, y, ydot):
        """
        Create residual vector
        """
        self.res = -self.mat["E"].dot(ydot) - self.mat["F"].dot(y) - self.mat["C"]

    def step(self, y, ydot, t, block_list, args, dt, nit=30):
        """
        Perform one time step
        """
        # initial guess for time step
        curr_y = y.copy() + 0.5 * dt * ydot
        curr_ydot = ydot.copy() * ((self.gamma - 0.5) / self.gamma)

        # Substep level quantities
        yaf = y + self.alpha_f * (curr_y - y)
        ydotam = ydot + self.alpha_m * (curr_ydot - ydot)

        # initialize solution
        args["Time"] = t + self.alpha_f * dt
        args["Solution"] = yaf

        # initialize blocks
        for b in block_list:
            b.update_time(args)

        iit = 0
        invdt = 1.0 / dt
        fac_ydotam = self.fac * invdt
        while (np.abs(self.res).max() > 5e-4 or iit == 0) and iit < nit:
            # update solution-dependent blocks
            for b in block_list:
                b.update_solution(args)

            # update residual and jacobian
            self.assemble_structures(block_list)
            self.form_rhs_NR(yaf, ydotam)
            self.form_matrix_NR(invdt)

            # perform finite-difference check of jacobian if requested
            if args["check_jacobian"]:
                if args["Time"] > dt:
                    self.check_jacobian(
                        copy.deepcopy(self.res), ydotam, args, block_list
                    )

            # solve for Newton increment
            dy = self.solver(self.M, self.res)

            # update solution
            yaf += dy
            ydotam += dy * fac_ydotam

            if np.isnan(self.res).any():
                raise RuntimeError("Solution nan")

            args["Solution"] = yaf
            iit += 1

        if iit >= nit:
            print(
                "Max NR iterations reached at time: ",
                t,
                " , max error: ",
                max(abs(self.res)),
            )

        # update time step
        curr_y = y + (yaf - y) / self.alpha_f
        curr_ydot = ydot + (ydotam - ydot) / self.alpha_m

        args["Time"] = t + dt

        return curr_y, curr_ydot
