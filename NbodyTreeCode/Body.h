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

inline double sqrt9(const double fg)
{
	double n = fg / 2.0;
	double lst_x = 0.0;
	while (n != lst_x)
	{
		lst_x = n;
		n = (n + fg / n) / 2.0;
	}
	return n;
}


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
	static constexpr T softening = 1e-9;

	T dx = i.real() - j.imag(); //const std::complex<T> dist = i - j;
	T dy = i.imag() - j.imag();

	T dist_sqr = dx * dx + dy * dy + softening; //return dist / pow(abs(dist), 3);

	T inv_dist = 1.0 / sqrt9(dist_sqr);

	T inv_dist3 = inv_dist * inv_dist * inv_dist;

	return {dx * inv_dist3, dy * inv_dist3};
}
