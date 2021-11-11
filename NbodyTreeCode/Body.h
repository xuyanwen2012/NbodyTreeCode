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

inline double sqrt9(const double x)
{
	union
	{
		int i;
		float x;
	} u;
	u.x = x;
	u.i = (1 << 29) + (u.i >> 1) - (1 << 22);

	// Two Babylonian Steps (simplified from:)
	// u.x = 0.5f * (u.x + x/u.x);
	// u.x = 0.5f * (u.x + x/u.x);
	u.x = u.x + x / u.x;
	u.x = 0.25f * u.x + x / u.x;

	return u.x;
}

inline double rsqrt64(const double number)
{
	const double x2 = number * 0.5;
	double y = number;
	uint64_t i = *reinterpret_cast<uint64_t*>(&y);
	i = 0x5fe6eb50c7b537a9 - (i >> 1);
	y = *reinterpret_cast<double*>(&i);
	y = y * (1.5 - (x2 * y * y));
	return y;
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

	T inv_dist = rsqrt64(dist_sqr); //T inv_dist = 1.0 / sqrt9(dist_sqr);

	T inv_dist3 = inv_dist * inv_dist * inv_dist;

	return {dx * inv_dist3, dy * inv_dist3};
}
