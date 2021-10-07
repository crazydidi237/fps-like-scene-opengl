#include "bridge.h"
using namespace vcl;

mesh initialize_plank() {
	mesh plank = mesh_load_file_obj("Assets/Village/Models/Props/OBJplank_01.obj");
	return plank;
}

mesh_drawable drawable_plank() {
	mesh_drawable plank=mesh_drawable(initialize_plank());
	plank.texture = opengl_texture_to_gpu(image_load_png("Assets/Village/Textures/Props/PNGprop_boardwalk_01_d.png"));
	plank.transform.scale = 0.4f;
	return plank;
}

mesh_drawable ground_plank() {
	mesh_drawable draw_plank = drawable_plank();
	draw_plank.transform.rotate = rotation({ 1,0,0 }, pi / 2);
	return draw_plank;
}

mesh_drawable transversal_ground_plank() {
	mesh_drawable draw_plank = drawable_plank();
	draw_plank.transform.rotate = rotation({ 0,0,1 }, pi / 2)*rotation({ 1,0,0 }, pi / 2);
	return draw_plank;
}

