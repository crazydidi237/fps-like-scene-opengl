#pragma once
#include "vcl/vcl.hpp"
#include "mesh_drawable_multitexture.h"


vcl::hierarchy_mesh_drawable create_fish_hierarchy(); 
std::vector<std::vector<vcl::vec3>> initialize_fishs(int number_of_fishs_per_line, float height, float x_start, float y_start);