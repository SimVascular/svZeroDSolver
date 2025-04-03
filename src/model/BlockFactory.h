// This file is part of svZeroDSolver licensed under Stanford University, The Regents of the University of 
//                                                   California, and others.
// 
// See the LICENSE.md file for license information
/**
 * @file BlockFactory.h
 * @brief Define the block factory functional
 */
#ifndef SVZERODSOLVER_MODEL_BLOCK_FACTORY_HPP_
#define SVZERODSOLVER_MODEL_BLOCK_FACTORY_HPP_

#include <functional>
#include <vector>

#include "Block.h"

/**
 * @brief General functional for the creation of different types of blocks
 */
using BlockFactoryFunc = std::function<Block*(int, Model*)>;

#endif
