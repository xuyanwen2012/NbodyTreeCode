#pragma once

#include <array>
#include <complex>

#ifdef __linux__
#include <algorithm>
#include <memory>
#endif

#include "Body.h"
#include "Rect.h"

namespace adaptive
{
	using vec2 = std::complex<double>;
	using body_ptr = std::shared_ptr<body<double>>;

	/// <summary>
	/// A node in the quadtree.
	/// </summary>
	struct tree_node
	{
		enum class direction { sw = 0, se, nw, ne };

		friend class quadtree;

		tree_node() = delete;

		tree_node(const int uid, const rect<double> bound, const size_t level)
			: uid(uid), level(level), bounding_box(bound), node_mass(0), is_leaf_(true)
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
		body_ptr content;

		/// <summary>
		///	 2 | 3
		/// ---+---
		///	 0 | 1
		/// </summary>
		std::array<tree_node*, 4> children{};

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
		bool is_empty() const { return content == nullptr; }

		/// <summary>
		///
		/// </summary>
		std::complex<double> center_of_mass() const { return weighted_pos / node_mass; }

	private:
		bool is_leaf_;
		void insert_body(const body_ptr& body_ptr);
		direction determine_quadrant(const vec2& pos) const;
		void split();
	};

	class quadtree
	{
	public:
		/// <summary>
		/// create a empty quadtree with only a node.
		/// </summary>
		quadtree();

		/// <summary>
		/// Use this function to insert a body into the quadtree.
		/// </summary>
		/// <param name="body_ptr"></param>
		void allocate_node_for_particle(const body_ptr& body_ptr);

		/// <summary>
		/// Once every particles are allocated into the quadtree, we can
		///	compute the center of masses and the quadtree is ready for
		///	inquiry.
		/// </summary>
		void compute_center_of_mass();

		/// <summary>
		/// 
		/// </summary>
		/// <param name="pos"></param>
		/// <returns></returns>
		std::complex<double> compute_force_at_iterative_dfs_array(const vec2& pos);

		// some statistical things
		size_t num_particles;
		static size_t num_nodes;
		static size_t depth;

	private:
		tree_node root_;

		bool check_theta(const tree_node* node, const vec2& pos) const;
		static inline std::complex<double> direct_compute(const body_ptr& body, const vec2& pos);
		static inline std::complex<double> estimate_compute(const tree_node* node, const vec2& pos);
	};
}
