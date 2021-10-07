#pragma once

#include "vcl/vcl.hpp"
#include"mesh_drawable_multitexture.h"

vcl::mesh initialize_second_terrain();
vcl::vec3 evaluate_mountain_pos(float u, float v);
vcl::vec3 evaluate_mountain_neg(float u, float v);
vcl::vec3 evaluate_main_terrain(float u, float v);
std::vector<vcl::mesh> create_main_terrain();
std::vector<vcl::mesh> create_mountain();

struct perlin_noise_parameters
{
	float persistency = 0.57f;
	float frequency_gain = 2.05f;
	int octave = 8;
	float terrain_height = 1.3f;
};

void update_terrain(vcl::mesh& mountain_pos, vcl::mesh_drawable& mountain_pos_visual,
    vcl::mesh& mountain_neg, vcl::mesh_drawable& mountain_neg_visual,
    vcl::mesh& main_terrain, vcl::mesh_drawable& main_terrain_visual,
	vcl::mesh& forest_terrain, vcl::mesh_drawable& forest_terrain_visual,
    vcl::mesh& path, vcl::mesh_drawable& path_visual,
    perlin_noise_parameters const& parameters, perlin_noise_parameters const& parameters_central, perlin_noise_parameters const& parameters_forest);
void update_rock(vcl::mesh& rock, vcl::mesh_drawable& rock_visual, perlin_noise_parameters const& parameters);
std::vector<vcl::vec3>  generate_houses_positions(int N, float safetyDistance, float safety_from_central_line);
void update_houses_positions(std::vector<vcl::vec3> &houses_positions, vcl::mesh& terrain);
vcl::mesh create_forest();
vcl::vec3 evaluate_forest_terrain(float u, float v);
std::vector<vcl::vec3>  generate_grass(int N, float safety_from_central_line);
void update_grass_positions(std::vector<vcl::vec3>& grass_positions, vcl::mesh& terrain);
std::vector<vcl::vec3>  generate_tree_positions(int N, float safetyDistance, float safety_from_central_line);
void update_tree_positions(std::vector<vcl::vec3>& tree_positions, vcl::mesh& terrain);
void update_camera_position(vcl::vec3& camera_position, vcl::mesh& terrain);