#include <iostream>
#include <complex>
#include <vector>
#include <random>

#include "AdaptiveQuadtree.h"
#include "Body.h"

#ifdef __linux__
#include <algorithm>
#include <cstdlib>
#endif

using namespace adaptive;

#ifdef _WIN32
double my_rand()
{
	static thread_local std::mt19937 generator; // NOLINT(cert-msc51-cpp)
	const std::uniform_real_distribution<double> distribution(0.0, 1.0);
	return distribution(generator);
}
#endif

#ifdef __linux__
double my_rand(const double f_min = 0.0, const double f_max = 1.0)
{
	const double f = static_cast<double>(rand()) / RAND_MAX;
	return f_min + f * (f_max - f_min);
}
#endif

/// <summary>
/// We have to extract this method out so that we can use the Princeton tool
/// on this function
/// </summary>
/// <param name="forces_n_log_n"></param>
/// <param name="qt"></param>
/// <param name="bodies"></param>
void estimate_forces(std::vector<vec2>& forces_n_log_n,
                     quadtree& qt,
                     const std::vector<body_ptr>& bodies)
{
	const auto num_bodies = bodies.size();
	for (size_t i = 0; i < num_bodies; ++i)
	{
		const auto force = qt.compute_force_at_iterative_dfs_array(bodies[i]->pos);
		forces_n_log_n.push_back(force);
	}
}

void _kernel_(quadtree& qt,
              const body_ptr& bodies)
{
	const auto result = qt.compute_force_at_iterative_dfs_array(bodies->pos);
	std::cout << result << std::endl;
}

/// <summary>
/// The Main Entry to my N-body program
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
/// <returns></returns>
int main(const int argc, char* argv[])
{
	static constexpr bool show_rmse = false;

	size_t num_bodies = 1024 * 1024;
	if (argc == 2)
	{
		num_bodies = atoi(argv[1]);
	}

	// The main particle table
	std::vector<body_ptr> bodies;

	// The force tables used to store results.
	std::vector<vec2> forces_n_squared;
	if (show_rmse)
	{
		forces_n_squared.reserve(num_bodies);
	}

	std::vector<vec2> forces_n_log_n;
	forces_n_log_n.reserve(num_bodies);

	// Initialization of positions/masses
	for (size_t i = 0; i < num_bodies; ++i)
	{
		const auto& pos = vec2{my_rand(), my_rand()};
		const auto& mass = my_rand() * 1.5;

		bodies.push_back(std::make_shared<body<double>>(i, pos, mass));
	}

	// -------- Do the N squared --------
	if (show_rmse)
	{
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
	}

	// -------- Do the NlogN --------
	auto qt = quadtree();

	// 1) Construct the Quadtree
	for (const auto& body : bodies)
	{
		qt.allocate_node_for_particle(body);
	}

	// 2) Calculate Centers of Mass
	qt.compute_center_of_mass();

	// 3) Estimate N-Body Forces
	estimate_forces(forces_n_log_n, qt, bodies);
	_kernel_(qt, bodies[0]);

	// -------- Do Analysis --------

	if (show_rmse)
	{
		vec2 tmp;
		for (size_t i = 0; i < num_bodies; ++i)
		{
			tmp += pow(forces_n_squared[i] - forces_n_log_n[i], 2);
		}

		const auto n = static_cast<double>(num_bodies);
		const auto rsme = sqrt(tmp / n);
		std::cout << "RSME = " << rsme << std::endl;
	}

	std::cout << "tree depth: " << quadtree::depth << std::endl;
	std::cout << "tree num nodes: " << quadtree::num_nodes << std::endl;

	return EXIT_SUCCESS;
}
