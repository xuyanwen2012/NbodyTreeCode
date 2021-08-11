#include <algorithm>
#include <array>
#include <iostream>
#include <random>
#include <complex>

#include "AdaptiveQuadtree.h"
#include "Body.h"

double my_rand()
{
	static thread_local std::mt19937 generator;
	const std::uniform_real_distribution distribution(0.0, 1.0);
	return distribution(generator);
}

int main()
{
	constexpr size_t num_bodies = 1024;

	// The main particle table
	std::vector<std::shared_ptr<body>> bodies;

	std::array<vec2, num_bodies> forces_n_squared;
	std::array<vec2, num_bodies> forces_n_log_n;

	// Initialization of positions/masses
	for (size_t i = 0; i < num_bodies; ++i)
	{
		const auto& pos = vec2{my_rand(), my_rand()};
		const auto& mass = my_rand() * 1.5;

		bodies.push_back(std::make_shared<body>(i, pos, mass));
	}

	// -------- Do the N squared --------
	for (size_t i = 0; i < num_bodies; ++i)
	{
		forces_n_squared[i] = {0, 0};
		for (size_t j = 0; j < num_bodies; ++j)
		{
			if (i == j)
			{
				continue;
			}

			const auto force = kernel_func(
				bodies[i]->pos,
				bodies[j]->pos
			);

			const auto fm = bodies[j]->mass * force;

			forces_n_squared[i] += fm;
		}
	}

	// -------- Do the NlogN --------
	auto qt = adaptive::quadtree();

	// 1) Construct the Quadtree
	for (const auto& body_ptr : bodies)
	{
		qt.allocate_node_for_particle(body_ptr);
	}

	// 2) Calculate Centers of Mass
	qt.compute_center_of_mass();

	// 3) Estimate N-Body Forces
	for (size_t i = 0; i < num_bodies; ++i)
	{
		forces_n_log_n[i] = qt.compute_force_at(bodies[i]->pos);
	}

	// -------- Do Analysis --------

	vec2 tmp;
	for (size_t i = 0; i < num_bodies; ++i)
	{
		tmp += pow(forces_n_squared[i] - forces_n_log_n[i], 2);
	}

	constexpr auto n = static_cast<double>(num_bodies);
	const auto rsme = sqrt(tmp / n);
	std::cout << "RSME = " << rsme << std::endl;

	return EXIT_SUCCESS;
}
