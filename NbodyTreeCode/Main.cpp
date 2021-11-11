#include <iostream>
#include <complex>
#include <vector>
#include <random>

#include "AdaptiveQuadtree.h"
#include "Body.h"

#ifdef __linux__
#include <algorithm>
#include <cstdlib>
#include <omp.h>

int DEC_NUM_CONSUMERS = 1;
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
	std::array<tree_node*, 1024> stack{};

	const auto num_bodies = bodies.size();
	for (size_t i = 0; i < num_bodies; ++i)
	{
		const auto force = qt.compute_force_at_iterative_dfs_array(stack, bodies[i]->pos, 1.0);
		forces_n_log_n[i] = - force;
	}
}


// ReSharper disable once CppInconsistentNaming
void _kernel_(std::vector<vec2>& forces_n_squared, // NOLINT(bugprone-reserved-identifier, clang-diagnostic-reserved-identifier)
              const std::vector<body_ptr>& bodies,
              const size_t num_bodies,
              const size_t num_to_sim,
              int t0,
              int t1)
{
	for (size_t i = 0; i < num_to_sim; ++i)
	{
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

/// <summary>
/// The Main Entry to my N-body program
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
/// <returns></returns>
int main(const int argc, char* argv[]) // NOLINT(bugprone-exception-escape)
{
	static constexpr bool show_rmse = false;

#ifdef __linux__
	srand(666);
#endif

	constexpr size_t num_bodies = 512;
	size_t num_to_sim = 512;

	if (argc == 2)
	{
		num_to_sim = std::stoi(argv[1]);
	}

	// The main particle table
	std::vector<body_ptr> bodies;

	// The force tables used to store results.
	std::vector<vec2> forces_n_squared(num_bodies);

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
		}
	}

	// 3) Estimate N-Body Forces
	_kernel_(forces_n_squared, bodies, num_bodies, num_to_sim, 0, 0);

	std::cout << "num_to_sim: " << num_to_sim << std::endl;

	return EXIT_SUCCESS;
}
