#include "AdaptiveQuadtree.h"

#include <iostream>
#include <list>
#include <algorithm>
#include <queue>

std::complex<double> adaptive::tree_node::get_gravity_at(const vec2& pos)
{
	std::complex<double> acc;

	if (is_leaf())
	{
		if (content->pos == pos) // making sure i != i
		{
			return 0;
		}

		// Direct computation
		const auto f = kernel_func(content->pos, pos);
		return content->mass * f;
	}

	const auto com = center_of_mass();
	const auto distance = com - pos;
	const auto norm = abs(distance);
	const auto geo_size = bounding_box.size.real();

	if (static double theta = 1.0; geo_size / norm < theta)
	{
		// we treat the quadtree cell as a source of long-range forces and use its center of mass.
		const auto f = kernel_func(com, pos);
		acc += node_mass * f;
	}
	else
	{
		// Otherwise, we will recursively visit the child cells in the quadtree.
		for (const auto child : children.value())
		{
			if (child->is_leaf() && child->is_empty())
			{
				continue;
			}

			acc += child->get_gravity_at(pos);
		}
	}

	return acc;
}

void adaptive::tree_node::insert_body(const std::shared_ptr<body>& body_ptr)
{
	if (is_leaf())
	{
		if (is_empty())
		{
			content = body_ptr;
			return;
		}

		// more than 1 particles are allocated into this node, need to split,
		// then re-insert the current content to the deeper levels
		split();

		const auto quadrant = static_cast<size_t>(determine_quadrant(content->pos));
		children->at(quadrant)->insert_body(content);

		content.reset();
	}

	const auto new_quadrant = static_cast<size_t>(determine_quadrant(body_ptr->pos));
	children->at(new_quadrant)->insert_body(body_ptr);
}

adaptive::tree_node::direction adaptive::tree_node::determine_quadrant(const vec2& pos) const
{
	const auto cx = bounding_box.center.real();
	const auto cy = bounding_box.center.imag();
	const auto x = pos.real();
	const auto y = pos.imag();

	if (x < cx)
	{
		if (y < cy)
		{
			return direction::sw;
		}
		return direction::nw;
	}
	if (y < cy)
	{
		return direction::se;
	}
	return direction::ne;
}

void adaptive::tree_node::split()
{
	const auto hw = bounding_box.size.real() / 2.0;
	const auto hh = bounding_box.size.imag() / 2.0;
	const auto cx = bounding_box.center.real();
	const auto cy = bounding_box.center.imag();

	const auto next_level = level + 1;

	const auto my_uid = uid * 10;

	const auto sw = new tree_node{my_uid + 0, rect{cx - hw / 2.0, cy - hh / 2.0, hw, hh}, next_level};
	const auto se = new tree_node{my_uid + 1, rect{cx + hw / 2.0, cy - hh / 2.0, hw, hh}, next_level};
	const auto nw = new tree_node{my_uid + 2, rect{cx - hw / 2.0, cy + hh / 2.0, hw, hh}, next_level};
	const auto ne = new tree_node{my_uid + 3, rect{cx + hw / 2.0, cy + hh / 2.0, hw, hh}, next_level};

	children = std::optional<std::array<tree_node*, 4>>{{sw, se, nw, ne}};
}

adaptive::quadtree::quadtree()
{
	num_particles_ = 0;
	root_ = tree_node(1, rect{0.5, 0.5, 1.0, 1.0}, 0);
}

void adaptive::quadtree::allocate_node_for_particle(const std::shared_ptr<body>& body_ptr)
{
	++num_particles_;
	root_.insert_body(body_ptr);
}

void adaptive::quadtree::compute_center_of_mass()
{
	std::queue<tree_node*> queue;
	std::vector<tree_node*> list;

	queue.push(&root_);
	while (!queue.empty())
	{
		const auto cur = queue.front();
		queue.pop();

		if (!cur->is_leaf())
		{
			for (auto child : cur->children.value())
			{
				queue.push(child);
			}
		}

		list.push_back(cur);
	}

	std::for_each(list.rbegin(), list.rend(),
	              [&](tree_node* node)
	              {
		              // sum the masses
		              double mass_sum = 0.0;
		              std::complex<double> weighted_pos_sum{0, 0};
		              if (node->is_leaf())
		              {
			              if (node->content != nullptr)
			              {
				              mass_sum = node->content->mass;
				              weighted_pos_sum = node->content->pos * node->content->mass;
			              }
		              }
		              else
		              {
			              for (const tree_node* child : node->children.value())
			              {
				              mass_sum += child->node_mass;
				              weighted_pos_sum += child->weighted_pos;
			              }
		              }

		              node->node_mass = mass_sum;
		              node->weighted_pos = weighted_pos_sum;
	              });
}

std::complex<double> adaptive::quadtree::compute_force_at(const vec2& pos)
{
	return root_.get_gravity_at(pos);
}
