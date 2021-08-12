#pragma once

#include <array>
#include <complex>
#include <optional>

#ifdef __linux__
#include <algorithm>
#include <memory>
#endif

#include "Body.h"
#include "Rect.h"

namespace adaptive
{
	using vec2 = std::complex<double>;

	/// <summary>
	/// A node in the quadtree.
	/// </summary>
	struct tree_node
	{
		enum class direction { sw = 0, se, nw, ne };

		friend class quadtree;

		tree_node() : uid(-1), level(0), node_mass(0)
		{
		}

		tree_node(const int uid, const rect<double> bound, const size_t level)
			: uid(uid), level(level), bounding_box(bound), node_mass(0)
		{
		}

		int uid;
		size_t level;

		/// <summary>
		/// I used center point as the position.
		///	Also the entire boarder of the whole universe is between [0..1]
		/// </summary>
		rect<double> bounding_box;

		/// <summary>
		///
		/// </summary>
		std::shared_ptr<body<double>> content;

		/// <summary>
		///	 2 | 3
		/// ---+---
		///	 0 | 1
		/// </summary>
		std::optional<std::array<tree_node*, 4>> children;

		/// <summary>
		/// This field stores the total mass of this node and its descendants
		/// </summary>
		double node_mass;

		/// <summary>
		///	The total sum of this node's and its descendants 'Position * mass'
		/// This is used to compute the center of mass, use it divide by 'node_mass'
		/// </summary>
		std::complex<double> weighted_pos;

		/// <summary>
		///
		/// </summary>
		[[nodiscard]] bool is_leaf() const { return !children.has_value(); }

		/// <summary>
		///
		/// </summary>
		[[nodiscard]] bool is_empty() const { return content == nullptr; }

		/// <summary>
		///
		/// </summary>
		[[nodiscard]] direction determine_quadrant(const vec2& pos) const;

		/// <summary>
		///
		/// </summary>
		[[nodiscard]] std::complex<double> center_of_mass() const { return weighted_pos / node_mass; }

		/// <summary>
		///
		/// </summary>
		std::complex<double> get_gravity_at(const vec2& pos);

	private:
		/// <summary>
		///
		/// </summary>
		///	<param name="body_ptr"> The body to be inserted into the quadtree. </param>
		void insert_body(const std::shared_ptr<body<double>>& body_ptr);

		/// <summary>
		///
		/// </summary>
		void split();
	};

	class quadtree
	{
	public:
		quadtree();

		void allocate_node_for_particle(const std::shared_ptr<body<double>>& body_ptr);
		void compute_center_of_mass();

		std::complex<double> compute_force_at_recursive(const vec2& pos);
		std::complex<double> compute_force_at_iterative_bfs(const vec2& pos);
		std::complex<double> compute_force_at_iterative_dfs(const vec2& pos);
		std::complex<double> compute_force_at_iterative_dfs_array(const vec2& pos);

		// some statistical things
		size_t num_particles;
		inline static size_t num_nodes = 1;
		inline static size_t depth = 0;

	private:
		tree_node root_;

		bool check_theta(const tree_node* node, const vec2& pos) const;
		static std::complex<double> direct_compute(const std::shared_ptr<body<double>>& body, const vec2& pos);
		static std::complex<double> estimate_compute(const tree_node* node, const vec2& pos);
	};
}
