/*
 * matrix.hpp
 *
 *  Created on: Jan 20, 2018
 *      Author: elchaschab
 */

#ifndef SRC_UTIL_MATRIX_HPP_
#define SRC_UTIL_MATRIX_HPP_

#include "geometry_util.hpp"
#include <array>
#include <vector>

typedef std::array<dim_t, 6> kernel_t;

class Matrix {
	std::vector<kernel_t> queue;
	kernel_t cache;
	bool cached = false;

public:
	kernel_t combine(const kernel_t& k1, const kernel_t& k2);
	bool isIdentity();
	Matrix& matrix(const kernel_t& m);
	Matrix& translate(const dim_t& tx, const dim_t& ty);
	Matrix& scale(const dim_t& sx, const dim_t& sy);
	Matrix& rotate(const dim_t& angle,const dim_t& rx,const dim_t& ry);
	Matrix& skewX(const dim_t& angle);
	Matrix& skewY(const dim_t& angle);
	kernel_t toArray();
	Point calc(const dim_t& x, const dim_t& y, bool isRelative);
};

#endif /* SRC_UTIL_MATRIX_HPP_ */
