#include "AdaptiveQuadtree.h"

#include <iostream>
#include <algorithm>
#include <queue>
#include <stack>

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

	static double theta = 1.0;
	if (geo_size / norm < theta)
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

void adaptive::tree_node::insert_body(const body_ptr& body_ptr)
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
	quadtree::depth = std::max(quadtree::depth, level);
	quadtree::num_nodes += 4;

	const auto my_uid = uid * 10;

	const auto sw = new tree_node{my_uid + 0, rect<double>{cx - hw / 2.0, cy - hh / 2.0, hw, hh}, next_level};
	const auto se = new tree_node{my_uid + 1, rect<double>{cx + hw / 2.0, cy - hh / 2.0, hw, hh}, next_level};
	const auto nw = new tree_node{my_uid + 2, rect<double>{cx - hw / 2.0, cy + hh / 2.0, hw, hh}, next_level};
	const auto ne = new tree_node{my_uid + 3, rect<double>{cx + hw / 2.0, cy + hh / 2.0, hw, hh}, next_level};

	children = tl::optional<std::array<tree_node*, 4>>{{sw, se, nw, ne}};
}

adaptive::quadtree::quadtree()
{
	num_particles = 0;
	root_ = tree_node(1, rect<double>{0.5, 0.5, 1.0, 1.0}, 0);
}

void adaptive::quadtree::allocate_node_for_particle(const body_ptr& body_ptr)
{
	++num_particles;
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
			for (auto child = cur->children.value().rbegin(); child != cur->children.value().rend(); ++child)
			{
			}

			std::for_each(cur->children.value().rbegin(),
			              cur->children.value().rend(),
			              [&](auto child)
			              {
			              });

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

std::complex<double> adaptive::quadtree::compute_force_at_recursive(const vec2& pos)
{
	return root_.get_gravity_at(pos);
}

std::complex<double> adaptive::quadtree::compute_force_at_iterative_bfs(const vec2& pos)
{
	std::complex<double> force;

	std::queue<tree_node*> queue;
	queue.push(&root_);

	while (!queue.empty())
	{
		std::cout<< queue.size() <<std::endl;

		const auto current = queue.front();
		queue.pop();

		if (current->is_leaf())
		{
			force += direct_compute(current->content, pos);
		}
		else if (check_theta(current, pos))
		{
			force += estimate_compute(current, pos);
		}
		else
		{
			// Otherwise, we will recursively visit the child cells in the quadtree.
			for (const auto child : current->children.value())
			{
				if (child->is_leaf() && child->is_empty()) // skip empty nodes
				{
					continue;
				}

				queue.push(child);
			}
		}
	}

	return force;
}

std::complex<double> adaptive::quadtree::compute_force_at_iterative_dfs(const vec2& pos)
{
	std::complex<double> force;

	std::stack<tree_node*> stack;
	stack.push(&root_);

	while (!stack.empty())
	{
		std::cout << stack.size() << std::endl;

		const auto current = stack.top();
		stack.pop();

		if (current->is_leaf())
		{
			force += direct_compute(current->content, pos);
		}
		else if (check_theta(current, pos))
		{
			force += estimate_compute(current, pos);
		}
		else
		{
			// Otherwise, we will recursively visit the child cells in the quadtree.
			for (const auto child : current->children.value())
			{
				if (child->is_leaf() && child->is_empty()) // skip empty nodes
				{
					continue;
				}

				stack.push(child);
			}
		}

	}

	return force;
}

std::complex<double> adaptive::quadtree::compute_force_at_iterative_dfs_array(const vec2& pos)
{
	std::complex<double> force;

	size_t stack_cp = 0;
	std::array<tree_node*, 1024> stack{};

	stack[++stack_cp] = &root_;

	while (stack_cp !=0)
	{
		const auto current = stack[stack_cp];
		stack[stack_cp--] = nullptr;

		if (current->is_leaf())
		{
			force += direct_compute(current->content, pos);
		}
		else if (check_theta(current, pos))
		{
			force += estimate_compute(current, pos);
		}
		else
		{
			// Otherwise, we will recursively visit the child cells in the quadtree.
			for (const auto child : current->children.value())
			{
				if (child->is_leaf() && child->is_empty()) // skip empty nodes
				{
					continue;
				}

				stack[++stack_cp] = child;
			}
		}

	}

	return force;
}

std::complex<double> adaptive::quadtree::direct_compute(const body_ptr& body, const vec2& pos)
{
	std::complex<double> force;

	if (body->pos != pos)
	{
		const auto f = kernel_func(body->pos, pos);
		force += f * body->mass;
	}

	return force;
}

bool adaptive::quadtree::check_theta(const tree_node* node, const vec2& pos) const
{
	const auto com = node->center_of_mass();
	const auto distance = com - pos;
	const auto norm = abs(distance);
	const auto geo_size = node->bounding_box.size.real();

	static double theta = 1.0;
	return geo_size / norm < theta;
}

std::complex<double> adaptive::quadtree::estimate_compute(const tree_node* node, const vec2& pos)
{
	const auto com = node->center_of_mass();
	const auto f = kernel_func(com, pos);
	return node->node_mass * f;
}

size_t adaptive::quadtree::depth = 0;

size_t adaptive::quadtree::num_nodes = 1;
