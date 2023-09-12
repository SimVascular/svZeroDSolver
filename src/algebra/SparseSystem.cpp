
#include "SparseSystem.h"
#include "../model/Model.h"

namespace algebra {

SparseSystem::SparseSystem() {}

SparseSystem::SparseSystem(unsigned int n) 
{
  F = Eigen::SparseMatrix<double>(n, n);
  E = Eigen::SparseMatrix<double>(n, n);
  D = Eigen::SparseMatrix<double>(n, n);
  C = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);

  jacobian = Eigen::SparseMatrix<double>(n, n);
  residual = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);
  dy = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(n);
}

SparseSystem::~SparseSystem() {}

void SparseSystem::clean() 
{
  // Cannot be in destructor because dynamically allocated pointers will be lost
  // when objects are assigned from temporary objects.
  //delete solver;
}

void SparseSystem::reserve(zd_model::Model *model) 
{
  auto num_triplets = model->get_num_triplets();
  F.reserve(num_triplets["F"]);
  E.reserve(num_triplets["E"]);
  D.reserve(num_triplets["D"]);
  model->update_constant(*this);
  model->update_time(*this, 0.0);

  Eigen::Matrix<double, Eigen::Dynamic, 1> dummy_y =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(residual.size());

  Eigen::Matrix<double, Eigen::Dynamic, 1> dummy_dy =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(residual.size());

  model->update_solution(*this, dummy_y, dummy_dy);

  F.makeCompressed();
  E.makeCompressed();
  D.makeCompressed();
  jacobian.reserve(num_triplets["F"] + num_triplets["E"]);  // Just an estimate
  update_jacobian(1.0);  // Update it once to have sparsity pattern
  jacobian.makeCompressed();
  solver->analyzePattern(jacobian);  // Let solver analyze pattern
}

void SparseSystem::update_residual( Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &ydot) 
{
  residual.setZero();
  residual -= C;
  residual.noalias() -= E * ydot;
  residual.noalias() -= F * y;
}

void SparseSystem::update_jacobian(double e_coeff) 
{
  jacobian.setZero();
  jacobian += F + D + E * e_coeff;
}

void SparseSystem::solve() 
{
  solver->factorize(jacobian);
  dy.setZero();
  dy += solver->solve(residual);
}

}; 

