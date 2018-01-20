#include "matrix.hpp"

kernel_t Matrix::combine(const kernel_t& k1, const kernel_t& k2) {
	return {
		k1[0] * k2[0] + k1[2] * k2[1],
		k1[1] * k2[0] + k1[3] * k2[1],
		k1[0] * k2[2] + k1[2] * k2[3],
		k1[1] * k2[2] + k1[3] * k2[3],
		k1[0] * k2[4] + k1[2] * k2[5] + k1[4],
		k1[1] * k2[4] + k1[3] * k2[5] + k1[5]
	};
}

bool Matrix::isIdentity() {
	if (!cached) {
		cache = toArray();
		cached = true;
	}

	kernel_t k = cache;

	if (k[0] == 1 && k[1] == 0 && k[2] == 0 && k[3] == 1 && k[4] == 0
			&& k[5] == 0)
		return true;
	else
		return false;
}

Matrix& Matrix::matrix(const kernel_t& m) {
	if (m[0] == 1 && m[1] == 0 && m[2] == 0 && m[3] == 1 && m[4] == 0
			&& m[5] == 0) {
		return *this;
	}
	cached = false;
	queue.push_back(m);
	return *this;
}

Matrix& Matrix::translate(const dim_t& tx, const dim_t& ty) {
	if (tx != 0 || ty != 0) {
		cached = false;
		queue.push_back( { 1, 0, 0, 1, tx, ty });
	}
	return *this;
}

Matrix& Matrix::scale(const dim_t& sx, const dim_t& sy) {
	if (sx != 1 || sy != 1) {
		cached = false;
		queue.push_back( { sx, 0, 0, sy, 0, 0 });
	}
	return *this;
}

Matrix& Matrix::rotate(const dim_t& angle, const dim_t& rx, const dim_t& ry) {
	dim_t rad, c, s;

	if (angle != 0) {
		translate(rx, ry);

		rad = angle * M_PI / 180;
		c = cos(rad);
		s = sin(rad);

		queue.push_back( { c, s, -s, c, 0, 0 });
		cached = false;

		translate(-rx, -ry);
	}
	return *this;
}

Matrix& Matrix::skewX(const dim_t& angle) {
	if (angle != 0) {
		cached = false;
		queue.push_back( { 1, 0, tan(angle * M_PI / 180), 1, 0, 0 });
	}
	return *this;
}

Matrix& Matrix::skewY(const dim_t& angle) {
	if (angle != 0) {
		cached = false;
		queue.push_back( { 1, tan(angle * M_PI / 180), 0, 1, 0, 0 });
	}
	return *this;
}

kernel_t Matrix::toArray() {
	if (cached) {
		return cache;
	}

	if (queue.empty()) {
		cache = {1, 0, 0, 1, 0, 0};
		return cache;
	}

	cache = queue[0];

	if (queue.size() == 1) {
		return cache;
	}

	for (size_t i = 1; i < queue.size(); i++) {
		cache = combine(cache, queue[i]);
	}

	return cache;
}

Point Matrix::calc(const dim_t& x, const dim_t& y, bool isRelative) {
	kernel_t m;
	// Don't change point on empty transforms queue
	if (queue.empty()) {
		return {x, y};
	}

	// Calculate final matrix, if not exists
	//
	// NB. if you deside to apply transforms to point one-by-one,
	// they should be taken in reverse order

	if (!cached) {
		cache = toArray();
	}

	m = cache;

	// Apply matrix to point
	return {
		x * m[0] + y * m[2] + (isRelative ? 0 : m[4]),
		x * m[1] + y * m[3] + (isRelative ? 0 : m[5])
	};
}

