#pragma once 

#include <cmath> 

inline bool is_equal(float x, float y)
{
	const float epsilon = 1.0e-10;
	return std::abs(x - y) <= epsilon * std::abs(x);
}

inline bool is_equal(double x, double y)
{
  const double epsilon = 1.0e-10;
  return std::abs(x - y) <= epsilon * std::abs(x);
}