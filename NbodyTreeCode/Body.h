#pragma once

#include <complex>

using vec2 = std::complex<double>;

struct body
{
	body(const int uid, const vec2 pos, const double mass)
		: uid(uid), pos(pos), mass(mass)
	{
	}

	int uid;
	vec2 pos;
	double mass;
};

inline vec2 kernel_func(const vec2& i, const vec2& j) { return log(abs(i - j)); }
