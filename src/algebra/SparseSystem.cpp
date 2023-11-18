// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "SparseSystem.h"

#include "Model.h"

SparseSystem::SparseSystem() {}

SparseSystem::SparseSystem(int n) {
  F = Eigen::SparseMatrix<double>(n, n);
  E = Eigen::SparseMatrix<double>(n, n);
  dC_dy = Eigen::SparseMatrix<double>(n, n);
  dC_dydot = Eigen::SparseMatrix<double>(n, n);
  C = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);

  jacobian = Eigen::SparseMatrix<double>(n, n);
  residual = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);
  dydot = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);
}

SparseSystem::~SparseSystem() {}

void SparseSystem::clean() {
  // Cannot be in destructor because dynamically allocated pointers will be lost
  // when objects are assigned from temporary objects.
  // delete solver;
}

void SparseSystem::reserve(Model *model) {
  auto num_triplets = model->get_num_triplets();
  F.reserve(num_triplets.F);
  E.reserve(num_triplets.E);
  dC_dy.reserve(num_triplets.D);
  dC_dydot.reserve(num_triplets.D);

  model->update_constant(*this);
  model->update_time(*this, 0.0);

  Eigen::Matrix<double, Eigen::Dynamic, 1> dummy_y =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(residual.size());

  Eigen::Matrix<double, Eigen::Dynamic, 1> dummy_dy =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(residual.size());

  model->update_solution(*this, dummy_y, dummy_dy);

  F.makeCompressed();
  E.makeCompressed();
  dC_dy.makeCompressed();
  dC_dydot.makeCompressed();
  jacobian.reserve(num_triplets.F + num_triplets.E);  // Just an estimate
  update_jacobian(1.0, 1.0);  // Update it once to have sparsity pattern
  jacobian.makeCompressed();
  solver->analyzePattern(jacobian);  // Let solver analyze pattern
}

void SparseSystem::update_residual(
    Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &ydot) {
  residual.setZero();
  residual -= C;
  residual.noalias() -= E * ydot;
  residual.noalias() -= F * y;
}

void SparseSystem::update_jacobian(double time_coeff_ydot,
                                   double time_coeff_y) {
  jacobian.setZero();
  jacobian += (E + dC_dydot) * time_coeff_ydot;
  jacobian += (F + dC_dy) * time_coeff_y;
}

void SparseSystem::solve() {
  solver->factorize(jacobian);
  dydot.setZero();
  dydot += solver->solve(residual);
}
