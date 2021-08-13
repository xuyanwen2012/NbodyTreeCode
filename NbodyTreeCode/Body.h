#pragma once

#include <complex>

template <typename T>
struct body
{
	body(const int uid, const std::complex<T> pos, const double mass)
		: uid(uid), pos(pos), mass(mass)
	{
	}

	int uid;
	std::complex<T> pos;
	double mass;
};

template <typename T>
std::complex<T> kernel_func(const std::complex<T>& i, const std::complex<T>& j)
{
	return log(abs(i - j));

	//const auto dist = (i - j) + 0.001;
	//return dist / ln(abs(dist));
}
