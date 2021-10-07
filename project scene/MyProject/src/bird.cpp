#include "bird.h"
#include "mesh_drawable_multitexture.h"
#include "terrain.hpp"

using namespace vcl;

hierarchy_mesh_drawable create_bird_hierarchy() {
	vec3 const ratio_body = vec3(0.21f, 0.12f, 0.12f);
	float const head_radius = 0.075f;
	float const length_arm = 0.15f;
	float const width_arm = 0.12f;

	hierarchy_mesh_drawable hierarchy;

	mesh_drawable body_bird = mesh_drawable(mesh_primitive_ellipsoid(ratio_body, vcl::vec3(0, 0, 0), 20, 20));
	mesh_drawable head_bird = mesh_drawable(mesh_primitive_sphere(head_radius, vec3(0, 0, 0), 20, 20));
	mesh_drawable eye = mesh_drawable(mesh_primitive_sphere(0.015f, { 0,0,0 }, 20, 20));
	eye.shading.color = { 0,0,0 };
	mesh_drawable beak = mesh_drawable(mesh_primitive_cone(0.045f, 0.075f, vec3(0, 0, 0), vec3(std::cos(3.14f / 8), 0, -std::sin(3.14f / 8)), false, 20, 20));
	beak.shading.color = { 0.3f, 0.2f, 0.7f };
	mesh_drawable left_arm = mesh_drawable(mesh_primitive_quadrangle(vec3(-width_arm, 0, 0), vec3(width_arm, 0, 0), vec3(0.5 * width_arm, length_arm, 0), vec3(-0.5 * width_arm, length_arm, 0)));
	mesh_drawable right_arm = mesh_drawable(mesh_primitive_quadrangle(vec3(-width_arm, 0, 0), vec3(width_arm, 0, 0), vec3(0.5 * width_arm, -length_arm, 0), vec3(-0.5 * width_arm, -length_arm, 0)));
	mesh_drawable left_fore_arm = mesh_drawable(mesh_primitive_quadrangle(vec3(-width_arm, 0, 0), vec3(width_arm, 0, 0), vec3(width_arm, length_arm, 0), vec3(-width_arm, length_arm, 0)));
	mesh_drawable right_fore_arm = mesh_drawable(mesh_primitive_quadrangle(vec3(-width_arm, 0, 0), vec3(width_arm, 0, 0), vec3(width_arm, -length_arm, 0), vec3(-width_arm, -length_arm, 0)));



	body_bird.texture = opengl_texture_to_gpu(image_load_png("Assets/Character/Textures/bird.png"), GL_REPEAT, GL_REPEAT);
	head_bird.texture = opengl_texture_to_gpu(image_load_png("Assets/Character/Textures/bird.png"), GL_REPEAT, GL_REPEAT);
	left_arm.texture = opengl_texture_to_gpu(image_load_png("Assets/Character/Textures/bird.png"), GL_REPEAT, GL_REPEAT);
	left_fore_arm.texture = opengl_texture_to_gpu(image_load_png("Assets/Character/Textures/bird.png"), GL_REPEAT, GL_REPEAT);
	right_fore_arm.texture = opengl_texture_to_gpu(image_load_png("Assets/Character/Textures/bird.png"), GL_REPEAT, GL_REPEAT);
	right_arm.texture = opengl_texture_to_gpu(image_load_png("Assets/Character/Textures/bird.png"), GL_REPEAT, GL_REPEAT);

	hierarchy.add(body_bird, "body");
	hierarchy.add(head_bird, "head", "body", { ratio_body[0] * 1.1, 0, ratio_body[2] * 0.4 });
	hierarchy.add(beak, "beak", "head", head_radius * vec3(1 / 2.0f, 0, -1 / 3.0f));
	hierarchy.add(eye, "left_eye", "head", head_radius * vec3(1 / 1.35f, 1 / 3.0f, 1 / 3.0f));
	hierarchy.add(eye, "right_eye", "head", head_radius * vec3(1 / 1.35f, -1 / 3.0f, 1 / 3.0f));
	hierarchy.add(left_fore_arm, "left_fore_arm", "body", { 0, ratio_body[1] - 0.03, 0 }); // extremity of the spherical body
	hierarchy.add(left_arm, "left_arm", "left_fore_arm", { 0, +length_arm,  0 });                        // the arm start at the center of the elbow
	hierarchy.add(right_fore_arm, "right_fore_arm", "body", { 0, -ratio_body[1] + 0.03, 0 }); // extremity of the spherical body
	hierarchy.add(right_arm, "right_arm", "right_fore_arm", { 0, -length_arm,  0 });                        // the arm start at the center of the elbow

	return hierarchy;
}

std::vector<std::vector<vcl::vec3>> initialize_birds(int number_of_birds_per_line, float height, float x_start, float y_start) {
	std::vector<vec3> positions, velocities, forces;
	positions.push_back({ x_start, y_start , height });
	velocities.push_back({ 0,0,0 });
	forces.push_back({ 0,0,0 });
	for (int i = 1; i < number_of_birds_per_line; i++) {
		positions.push_back({x_start-i, y_start + i, height});
		positions.push_back({ x_start -i, y_start  -i, height});
		velocities.push_back({ 0,0,0 });
		velocities.push_back({ 0,0,0 });
		forces.push_back({ 0,0,0 });
		forces.push_back({ 0,0,0 });
	}

	return { positions, velocities, forces };
}