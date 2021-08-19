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

/// <summary>
/// The main kernel function used to compute the pairwise force between
///	particles. 
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="i">The first particle.</param>
/// <param name="j">The second particle.</param>
/// <returns>The force.</returns>
template <typename T>
std::complex<T> kernel_func(const std::complex<T>& i, const std::complex<T>& j)
{
	const auto dist = i - j;
	return dist / pow(abs(dist), 3);
}
