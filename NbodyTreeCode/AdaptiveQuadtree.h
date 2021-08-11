#pragma once

#include <array>
#include <array>
#include <complex>
#include <optional>

#include "Body.h"
#include "Rect.h"

namespace adaptive
{
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
		std::shared_ptr<body> content;

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
		_NODISCARD bool is_leaf() const { return !children.has_value(); }

		/// <summary>
		///
		/// </summary>
		_NODISCARD bool is_empty() const { return content == nullptr; }

		/// <summary>
		///
		/// </summary>
		_NODISCARD direction determine_quadrant(const vec2& pos) const;

		/// <summary>
		///
		/// </summary>
		_NODISCARD std::complex<double> center_of_mass() const { return weighted_pos / node_mass; }

		/// <summary>
		///
		/// </summary>
		std::complex<double> get_gravity_at(const vec2& pos);
	private:
		/// <summary>
		///
		/// </summary>
		///	<param name="body_ptr"> The body to be inserted into the quadtree. </param>
		void insert_body(const std::shared_ptr<body>& body_ptr);

		/// <summary>
		///
		/// </summary>
		void split();
	};

	class quadtree
	{
	public:
		quadtree();
		void allocate_node_for_particle(const std::shared_ptr<body>& body_ptr);
		void compute_center_of_mass();
		std::complex<double> compute_force_at(const vec2& pos);

	private:
		tree_node root_;
		size_t num_particles_;
	};
}
