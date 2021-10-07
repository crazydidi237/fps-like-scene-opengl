#pragma once

#include "vcl/vcl.hpp"


vcl::mesh create_tree_trunk_cylinder(float r, float h, vcl::vec3 Deviation, vcl::vec3 Pos_init);
vcl::mesh create_leaf(vcl::vec3 Deviation, vcl::vec3 Pos_init);
void tree_rec(float r, float h, vcl::vec3 Deviation, vcl::vec3 Pos_init, vcl::mesh& Abr, int n);
std::vector<vcl::vec3> Foliage_pos();
vcl::mesh create_tree();