#pragma once

#include <complex>

template <typename T>
struct rect
{
	rect() = default;

	rect(const T cx, const T cy, const T w, const T h) : center(cx, cy), size(w, h)
	{
	}

	std::complex<T> center;
	std::complex<T> size;
};
