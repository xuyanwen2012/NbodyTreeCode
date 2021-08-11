#pragma once

#include <complex>

template <typename T,
          typename = std::enable_if_t<std::is_floating_point_v<T>, T>>
struct rect
{
	rect() = default;

	rect(const T cx, const T cy, const T w, const T h) : center(cx, cy), size(w, h)
	{
	}

	std::complex<T> center;
	std::complex<T> size;
};
