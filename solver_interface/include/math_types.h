/* This file is originally part of flightmare and under MIT License

* Original source: https://github.com/uzh-rpg/flightmare/blob/master/flightlib/include/flightlib/common/types.hpp
* Copyright (c) 2021 Robotics and Perception Group, 
* University of Zurich, Switzerland
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.

*/

#pragma once

#include <Eigen/Dense>
#include <Eigen/Core>
#include "../extern/solver/PerchingSolver/include/PerchingSolver.h"

// Define the scalar type used.
using Scalar = PerchingSolver_float;
static constexpr Scalar INF = std::numeric_limits<Scalar>::infinity();

// Define `Dynamic` matrix size.
static constexpr int Dynamic = Eigen::Dynamic;

// Using shorthand for `Matrix<rows, cols>` with scalar type.
template<int rows = Dynamic, int cols = rows>
using Matrix = Eigen::Matrix<Scalar, rows, cols>;

// Using shorthand for `Vector<rows>` with scalar type.
template<int rows = Dynamic>
using Vector = Matrix<rows, 1>;

// Using `Quaternion` with type.
using Quaternion = Eigen::Quaternion<Scalar>;

// Using `Ref` for modifier references.
template<class Derived>
using Ref = Eigen::Ref<Derived>;

// Using `ConstRef` for constant references.
template<class Derived>
using ConstRef = const Eigen::Ref<const Derived>;