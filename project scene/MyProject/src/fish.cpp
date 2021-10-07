#include "fish.h"
#include "terrain.hpp"
#include "mesh_drawable_multitexture.h"

using namespace vcl;

hierarchy_mesh_drawable create_fish_hierarchy() {
	hierarchy_mesh_drawable hierarchy;
	// The geometry of the body is an ellipsoid
	mesh_drawable body = mesh_drawable(mesh_primitive_ellipsoid({ 0.1,0.5,0.25 }));
	// Geometry of the eyes: black spheres
	mesh_drawable eye = mesh_drawable(mesh_primitive_sphere(0.05f, { 0,0,0 }, 20, 20));
	eye.shading.color = { 0,0,0 };

	//tail as a triangle
	mesh_drawable tail = mesh_drawable(mesh_primitive_triangle({ 0,-0.3,-0.3 }, { 0,0,0 }, { 0,-0.3,0.3 }));

	// fins displayed as triangles
	mesh_drawable shoulder_left = mesh_drawable(mesh_primitive_triangle({ -0.05,0.05,0 }, { 0,0,0 }, { -0.05,-0.05,0 }));

	mesh_drawable shoulder_right = mesh_drawable(mesh_primitive_triangle({ 0.05,0.05,0 }, { 0,0,0 }, { 0.05,-0.05,0 }));

	// Build the hierarchy:
	// ------------------------------------------- //
	// Syntax to add element
	//   hierarchy.add(visual_element, element_name, parent_name, (opt)[translation, rotation])
	// The root of the hierarchy is the body
	hierarchy.add(body, "body");
	// Eyes positions are set 
	hierarchy.add(eye, "eye_left", "body", 0.25 * vec3(1 / 5.0f, 1.5f, 0.3f));
	hierarchy.add(eye, "eye_right", "body", 0.25 * vec3(-1 / 5.0f, 1.5f, 0.3f));
	// Set the left part 
	hierarchy.add(shoulder_left, "shoulder_left", "body", { -0.1f,0.1,0 }); // extremity of the spherical body
	hierarchy.add(tail, "tail", "body", { 0.f,-0.5,0. });
	// Set the right part
	hierarchy.add(shoulder_right, "shoulder_right", "body", { 0.1f,0.1,0 });                    

	return hierarchy;
}

std::vector<std::vector<vcl::vec3>> initialize_fishs(int number_of_fishs_per_line, float height, float x_start, float y_start) {

	std::vector<vec3> positions, velocities, forces;
	for (int i = 0; i < number_of_fishs_per_line; i++) {
		for (int j = 0; j < 3*number_of_fishs_per_line; j++) {
			positions.push_back({ x_start + i, y_start +  1.5*j, height });
			velocities.push_back({ 0,0,0 });
			forces.push_back({ 0,0,0 });
		}	
	}

	return { positions, velocities, forces };
}