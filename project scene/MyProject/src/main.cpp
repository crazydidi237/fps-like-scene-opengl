#include "vcl/vcl.hpp"
#include <iostream>
#include <list>

#include "terrain.hpp"
#include "bridge.h"
#include "mesh_drawable_multitexture.h"
#include "mesh_drawable_multitexture2.h"
#include "tree.h"
#include "bird.h"
#include "environment_map.hpp"
#include "water_fbo.hpp"
#include "mesh_drawable_multitexture3.hpp"
#include <fish.h>

//#include "tree.hpp"

using namespace vcl;

struct gui_parameters {
	bool display_frame = true;
	bool add_sphere = true;
	bool display_wireframe = false;
};


struct keyboard_state_parameters {
	bool left = false;
	bool right = false;
	bool up = false;
	bool down = false;
};

struct user_interaction_parameters {
	vec2 mouse_prev;
	timer_fps fps_record;
	mesh_drawable global_frame;
	gui_parameters gui;
	bool cursor_on_gui;
	keyboard_state_parameters keyboard_state;
	bool mouse_click;
};
user_interaction_parameters user;

struct scene_environment
{
	camera_around_center camera;
	mat4 projection;
	vec3 light;
	float t;
	vec4 clip_plan;
	vec4 clip_plan1;
	vec4 clip_plan2;//clipping planes for water reflection and refraction texture
	vec3 cam_Pos; //Uniform cam position
};
scene_environment scene;
scene_environment scene_fps;
scene_environment *scene_to_draw(&scene);


// water variables
mesh_drawable_multitexture3 water_plane;
float Movefact = 0;
float const Wavespeed = 0.03f;
timer_basic timer1;

mesh_drawable ground;

//camera variables
vec3 cameraPosition = evaluate_main_terrain(0.5, 0) + vec3(0,0,1);
vec3 cameraFront = vec3(0.0f, 1.0f, 0.0f);
vec3 cameraUp = vec3(0.0f, 0.0f, 1.0f);

float yaw = 0.0f;
float pitch = 0.0f;
float fov_cam = 50.0f * pi / 180.0f;

bool first_mouse = true;
float lastX = 1280.0f / 2.0;
float lastY = 400.0f / 2.0;

bool camera_fps_active = false;
//end camera variables

perlin_noise_parameters parameters;
perlin_noise_parameters parameters_central;
perlin_noise_parameters parameters_forest;

struct particle_structure
{
	vcl::vec3 p; // Position
	vcl::vec3 v; // Speed
};

void keyboard_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void mouse_move_callback(GLFWwindow* window, double xpos, double ypos);
void window_size_callback(GLFWwindow* window, int width, int height);

void initialize_data();
void display_frame();
void display_interface();
void update_camera();

mesh_drawable cube;

mesh_drawable_multitexture mult_test;
mesh_drawable_multitexture main_terrain;
mesh_drawable_multitexture mountain_pos;
mesh_drawable_multitexture mountain_neg;
mesh_drawable_multitexture forest_terrain;
mesh_drawable forest_pos;
mesh_drawable forest_neg;
mesh_drawable tree;
mesh_drawable mushroom;
std::vector<mesh_drawable> houses;
std::vector<mesh_drawable> rock;
mesh_drawable path_in_snow;
mesh_drawable plank; //sinlge plank for the construction of the bridge
mesh_drawable single_ground_plank;
mesh_drawable single_ground_transversal_plank;
mesh_drawable single_ground_longitudinal_plank;

mesh_drawable billboard_fire;
mesh mesh_billboard_fire;

mesh_drawable tree1;

mesh mesh_main_terrain;
mesh mesh_mountain_pos;
mesh mesh_mountain_neg;
mesh mesh_forest;
mesh mesh_path_in_snow;


std::vector<mesh> mesh_houses;
std::vector<mesh> mesh_vegetation;
std::vector<mesh> mesh_rock;
mesh_drawable skybox_drawable;

std::vector<GLuint> house_textures;
std::vector<GLuint> rock_textures;
GLuint vegetation_textures;


//mesh rock;
mesh_drawable rock_drawable;
std::vector<mesh_drawable> textured_tree;

std::vector<vcl::vec3> tree_position;
std::vector<vcl::vec3> mushroom_position;

std::vector<vcl::vec3> houses_positions;
std::vector<vcl::vec3> vegetation_position;
std::vector<vcl::vec3> rock_position;

std::vector<int> house_indices;
std::vector<int> vegetation_indices;
std::vector<int> rock_indices;

std::vector<mat3> houses_rand_rotation;
std::vector<float> houses_rand_scale;

std::vector<mat3> vegetation_rand_rotation;
std::vector<float> vegetation_rand_scale;

std::vector<mat3> rock_rand_rotation;
std::vector<float> rock_rand_scale;

std::vector<mat3> tree_rand_rotation;
std::vector<float> tree_rand_scale;

mat3 initial_house_rot;

float height_bridge_handle = 0.8f;
float x_base_bridge = 14.2f, z_base_bridge = 1.4f;

int number_of_houses = 50;
int number_of_vegetation = 120;
int number_of_rock = 60;
int number_of_grass = 50;

std::list<particle_structure> snow_particles; // Storage of all currently active particles
mesh_drawable snow_sphere;
timer_event_periodic timer(0.08f);

int i = 0, j = 0;

mesh test;

vec3 memory1;
vec3 memory2;
vec3 memory3;
vec3 memory4;

mesh mesh_billboard;
std::vector<vec3> billboard_positions ;
mesh_drawable billboard;
float t;
vec3 wind;
vec3 wind_flag;
float friction;
float last_timer;
timer_event_periodic timer_billboard(10.0f);

timer_event_periodic timer_enemy(10.0f);



mesh mesh_flag;
mesh mesh_flag_support;
mesh_drawable flag_support;
mesh_drawable flag;
mesh_drawable flag2;
int N_flag;
std::vector<std::vector<vcl::vec3>> velocities;
float L0;
float L0_array;
float L0_cisaillement;
float L0_flexion;


mesh mesh_billboard_grass;
std::vector<vec3> billboard_positions_grass;
std::vector<float> billboard_scale_grass;
mesh_drawable billboard_grass;


std::vector<vcl::vec3> velocities_animals;
std::vector<vcl::vec3> positions_animals;
std::vector<vcl::vec3> positions_to_avoid;
float r0 = 0.2f;
float V = 0.01f;
float r_bird = 0.5f;
float V_bird = 0.1f;
float m_bird = 0.05f;

float r_fish = 0.5f;
float V_fish = 0.1f;
float m_fish = 0.05f;

std::vector<vec3> destination_animals;
std::vector<bool> new_destination;


hierarchy_mesh_drawable bird_hierarchy;
std::vector<vec3> bird_positions;
std::vector<vec3> bird_velocities;
std::vector<vec3> bird_forces;

hierarchy_mesh_drawable fish_hierarchy;
std::vector<vec3> fish_positions;
std::vector<vec3> fish_velocities;
std::vector<vec3> fish_forces;

hierarchy_mesh_drawable enemy_hierarchy;
std::vector<vec3> enemy_positions;
std::vector<vec3> enemy_velocities;
std::vector<vec3> enemy_forces;

mesh_drawable mouse_pointer;
bool gun_animation;
float t_clicked;


mesh_drawable sphere;

mesh_drawable gun;


void init_camera_fps(mesh& terrain) {
	update_camera_position(cameraPosition, terrain);
	scene_fps.camera.look_at(cameraPosition, cameraPosition + cameraFront, cameraUp);
}


vec3 spring_force(vec3& p_i, vec3& p_j, float L_0, float K)
{
	vec3 p = p_i - p_j;
	float L = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);

	// TO DO: correct the computation of this force value
	return -K * (L - L_0) * p / L;
}

vec3 lennard_jones_forces(vec3& p_i, std::vector<vec3>& positions, float r0, float V) {

	vec3 force = vec3(0, 0, 0);
	
		for (auto p : positions) {
			vec3 ptemp = p - p_i;
			float r = norm(ptemp);

			if (r != 0) {
				float u = V * (pow((r0 / r), 12) - 2 * pow((r0 / r), 6));
				force += -u * ptemp / (r * r);

			}


		}

	
	return force;

}

std::vector<std::vector<vcl::vec3>> enemy_spawner(std::vector<vcl::vec3> positions, std::vector<vcl::vec3> velocities, std::vector<vcl::vec3> forces, std::vector<vcl::vec3> positions_spawner, int rate) {
	int N = positions_spawner.size();
	
	for (int i = 0; i < rate; i++) {
		positions.push_back(positions_spawner[rand_interval(0, N - 1)] +vec3(0,0,0.3f));
		velocities.push_back({ 0,0,0 });
		forces.push_back({ 0,0,0 });
	}
	std::cout << "spawning\n";

	return { positions, velocities, forces };
}

std::vector<std::vector<vcl::vec3>> animate_enemy(hierarchy_mesh_drawable& hierarchy, std::vector<vcl::vec3> positions, std::vector<vcl::vec3> velocities, std::vector<vcl::vec3> forces, float m, float mu, float dtf, float r_fish, float V_fish, float freq) {
	//Rotation of the head
	hierarchy["head"].transform.rotate = rotation({ 0,1,0 }, 0.2 * std::sin(10 * 3.14f * (t - 0.4f)));
	// Rotation of the shoulder-left around the y axis
	hierarchy["left_fore_arm"].transform.rotate = rotation({ -1,0,0 }, std::sin(10 * 3.14f * (t - 0.4f)));
	// Rotation of the arm-left around the y axis (delayed with respect to the shoulder)
	hierarchy["left_arm"].transform.rotate = rotation({ -1,0,0 }, std::sin(10 * 3.14f * (t - 0.6f)));

	// Rotation of the shoulder-right around the y axis
	hierarchy["right_fore_arm"].transform.rotate = rotation({ 1,0,0 }, std::sin(10 * 3.14f * (t - 0.4f)));
	// Rotation of the arm-right around the y axis (delayed with respect to the shoulder)
	hierarchy["right_arm"].transform.rotate = rotation({ 1,0,0 }, std::sin(10 * 3.14f * (t - 0.6f)));

	for (int i = 0; i < forces.size(); i++) {
		forces[i] = lennard_jones_forces(positions[i], positions, r_bird, V_bird);
	}

	
	for (int i = 0; i < velocities.size(); i++) {
		vec3 p_i = positions[i];
		velocities[i] = (1 - mu) * velocities[i] + 0.02 * dtf * forces[i] / (m)+vec3(0.0, 0, 0.5 * std::sin(freq * t))+0.7*t*normalize(scene_to_draw->camera.position() - positions[i]);
		positions[i] = positions[i] + 0.02 * dtf * velocities[i];
		hierarchy["body"].transform.translate = vec3(0, 0, 0.05f * (1 + std::sin(2 * 3.14f * t))) + positions[i];
		hierarchy["body"].transform.rotate = scene_to_draw->camera.orientation_camera * rotation({ 0,0,1 }, -pi / 2) * rotation({ 0,1,0 }, -pi / 2);
		//vec3 c = normalize( cross({ 1,0,0 }, positions[i] - p_i));
		//hierarchy["body"].transform.rotate = rotation(c, std::acos(dot({ 1,0,0 }, positions[i] - p_i)/norm(positions[i] - p_i)) );
		hierarchy.update_local_to_global_coordinates();
		draw(hierarchy, *scene_to_draw);
	}

	return { positions, velocities, forces };

}


void create_ground_bridge(float x_offset, scene_environment& scene) {
	for (int i = 0; i < 13; i++) {
		single_ground_plank.transform.translate = { x_base_bridge + 0.14 * i + x_offset, 0, z_base_bridge };
		draw(single_ground_plank, scene);
	}//0.35f pour l'incrément et (2.35f, -2.35f) pour les bois de coté sans le facteur scale
	single_ground_transversal_plank.transform.translate = { x_base_bridge + 0.14 * 6 + x_offset, 0.94f, z_base_bridge };
	draw(single_ground_transversal_plank, scene);
	single_ground_transversal_plank.transform.translate = { x_base_bridge + 0.14 * 6 + x_offset, -0.94f, z_base_bridge };
	draw(single_ground_transversal_plank, scene);
	single_ground_transversal_plank.transform.translate = { x_base_bridge + 0.14 * 6 + x_offset, -0.94f, z_base_bridge + height_bridge_handle };
	draw(single_ground_transversal_plank, scene);
	single_ground_transversal_plank.transform.translate = { x_base_bridge + 0.14 * 6 + x_offset, 0.94f, z_base_bridge + height_bridge_handle };
	draw(single_ground_transversal_plank, scene);

	single_ground_longitudinal_plank.transform.translate = { x_base_bridge + 0.14 * 13 + x_offset - 0.04f, 0.94f, z_base_bridge };
	draw(single_ground_longitudinal_plank, scene);
	single_ground_longitudinal_plank.transform.translate = { x_base_bridge + 0.14 * 13 + x_offset - 0.04f, -0.94f, z_base_bridge };
	draw(single_ground_longitudinal_plank, scene);

}

std::vector<std::vector<vcl::vec3>> animate_bird_hierarchy(hierarchy_mesh_drawable& hierarchy, std::vector<vcl::vec3> positions, std::vector<vcl::vec3> velocities, std::vector<vcl::vec3> forces, float m, float mu, float dtf,float r_bird, float V_bird, float freq) {
	
	//Rotation of the head
	hierarchy["head"].transform.rotate = rotation({ 0,1,0 }, 0.2 * std::sin(10 * 3.14f * (t - 0.4f)));
	// Rotation of the shoulder-left around the y axis
	hierarchy["left_fore_arm"].transform.rotate = rotation({ -1,0,0 }, std::sin(10 * 3.14f * (t - 0.4f)));
	// Rotation of the arm-left around the y axis (delayed with respect to the shoulder)
	hierarchy["left_arm"].transform.rotate = rotation({ -1,0,0 }, std::sin(10 * 3.14f * (t - 0.6f)));

	// Rotation of the shoulder-right around the y axis
	hierarchy["right_fore_arm"].transform.rotate = rotation({ 1,0,0 }, std::sin(10 * 3.14f * (t - 0.4f)));
	// Rotation of the arm-right around the y axis (delayed with respect to the shoulder)
	hierarchy["right_arm"].transform.rotate = rotation({ 1,0,0 }, std::sin(10 * 3.14f * (t - 0.6f)));

	for (int i = 0; i < forces.size(); i++) {
		forces[i] = lennard_jones_forces(positions[i], positions, r_bird, V_bird);
	}

	if (bird_positions[8][0] >= 50) {

		positions.clear();
		velocities.clear();
		forces.clear();
		std::vector<std::vector<vec3>> temp_birds = initialize_birds(5, 10, -25, -3);
		std::vector<std::vector<vec3>> temp_birds1 = initialize_birds(4, 8, -15.0, 5.0f);
		positions = temp_birds[0];
		positions.insert(positions.end(), temp_birds1[0].begin(), temp_birds1[0].end());
		velocities = temp_birds[1];
		velocities.insert(velocities.end(), temp_birds1[1].begin(), temp_birds1[1].end());
		forces = temp_birds1[2];
		forces.insert(forces.end(), temp_birds1[2].begin(), temp_birds1[2].end());
		return { positions, velocities, forces };
	}

	for (int i = 0; i < velocities.size(); i++) {
		vec3 p_i = positions[i];
		float angle;
		velocities[i] = (1 -  mu) * velocities[i] + 0.02*dtf * forces[i] / (m) + vec3(0.7*t, 0, 0.5*std::sin(freq*t)) ;
		positions[i] = positions[i] + 0.02*dtf * velocities[i];
		hierarchy["body"].transform.translate = vec3( 0,0,0.5f * (1 + std::sin(2 * 3.14f * t)) ) + positions[i];
		//vec3 c = normalize( cross({ 1,0,0 }, positions[i] - p_i));
		//hierarchy["body"].transform.rotate = rotation(c, std::acos(dot({ 1,0,0 }, positions[i] - p_i)/norm(positions[i] - p_i)) );
		hierarchy.update_local_to_global_coordinates();
		draw(hierarchy, *scene_to_draw);
	}

	return { positions, velocities, forces };


	
}

std::vector<std::vector<vcl::vec3>> animate_fish_hierarchy(hierarchy_mesh_drawable& hierarchy, std::vector<vcl::vec3> positions, std::vector<vcl::vec3> velocities, std::vector<vcl::vec3> forces, float m, float mu, float dtf, float r_fish, float V_fish, float freq) {

	
	// Update the current time
	timer.update();
	float const t = timer.t;

	
	// Rotation of the fin-left around the y axis
	hierarchy["shoulder_left"].transform.rotate = rotation({ 0,1,0 }, std::sin(2 * 3.14f * (t - 0.4f)));
	// Rotation of the fin-right around the y axis
	hierarchy["shoulder_right"].transform.rotate = rotation({ 0,-1,0 }, std::sin(2 * 3.14f * (t - 0.4f)));
	//tail rotation
	hierarchy["tail"].transform.rotate = rotation({ 0,0,1 }, std::sin(2 * 3.14f * (t - 0.4f)));
	// update the global coordinates
	hierarchy.update_local_to_global_coordinates();

	// display the hierarchy
	
	for (int i = 0; i < forces.size(); i++) {
		forces[i] = lennard_jones_forces(positions[i], positions, r_fish, V_fish);
	}

	if (positions[8*2][1] >= water_plane.transform.translate.y + 30) {

		positions.clear();
		velocities.clear();
		forces.clear();
		std::vector<std::vector<vec3>> temp_fishs = initialize_fishs(3, ground.transform.translate.z + 0.1, water_plane.transform.translate.x - 2, water_plane.transform.translate.y - 30);
		std::vector<std::vector<vec3>> temp_fishs1 = initialize_fishs(4, ground.transform.translate.z + 0.1, water_plane.transform.translate.x + 3, water_plane.transform.translate.y - 25);
		positions = temp_fishs[0];
		positions.insert(positions.end(), temp_fishs1[0].begin(), temp_fishs1[0].end());
		velocities = temp_fishs[1];
		velocities.insert(velocities.end(), temp_fishs1[1].begin(), temp_fishs1[1].end());
		forces = temp_fishs[2];
		forces.insert(forces.end(), temp_fishs1[2].begin(), temp_fishs1[2].end());
		return { positions, velocities, forces };
	}

	for (int i = 0; i < velocities.size(); i++) {
		vec3 p_i = positions[i];
		float angle;
		velocities[i] = (1 - mu) * velocities[i] + 0.02*dtf * forces[i] / (m)+vec3(0, 0.1 * t, 0.5 * std::sin(freq * t));
		positions[i] = positions[i] + 0.02*dtf * velocities[i];
		hierarchy["body"].transform.translate = vec3(0.3f * (1 + std::sin(2 * 3.14f * t)), 0, 0) + positions[i];
		hierarchy.update_local_to_global_coordinates();
		draw(hierarchy, *scene_to_draw);
		//vec3 c = normalize( cross({ 1,0,0 }, positions[i] - p_i));
		//hierarchy["body"].transform.rotate = rotation(c, std::acos(dot({ 1,0,0 }, positions[i] - p_i)/norm(positions[i] - p_i)) );
		
	}

	return { positions, velocities, forces };
	

}

void create_full_ground_bridge(int N, scene_environment& scene) {
	for (int i = 0; i < N; i++) {
		create_ground_bridge(i * (0.14 * 13), scene);
	}

	single_ground_longitudinal_plank.transform.translate = { x_base_bridge , 0.94f, z_base_bridge };
	draw(single_ground_longitudinal_plank, scene);
	single_ground_longitudinal_plank.transform.translate = { x_base_bridge , -0.94f, z_base_bridge };
	draw(single_ground_longitudinal_plank, scene);
}


void initialize_flag(vec3 position_flag) {
	mesh_flag_support = mesh_primitive_cylinder(0.02f, { 0,0,0 }, { 0,0,3 });
	mesh_flag_support.push_back(mesh_primitive_cylinder(0.02f, { 0,0,3 }, { 0,1.0,3 }));
	flag_support = mesh_drawable(mesh_flag_support);


	mesh_flag = mesh_primitive_grid();
	flag = mesh_drawable(mesh_flag);
	flag.transform.translate = position_flag;
	N_flag = sqrt(mesh_flag.position.size());
	for (int i = 0; i < N_flag; i++) {
		std::vector<vec3> temp2;
		velocities.push_back(temp2);
		for (int j = 0; j < N_flag; j++) {
			velocities[i].push_back({ 0,0,0 });
		}
	}

	L0 = 0.1f; // Rest length between A and B
	L0_cisaillement = L0 * sqrt(2);
	L0_flexion = 2 * L0;

	flag.texture = opengl_texture_to_gpu(image_load_png("Assets/Village/Textures/Flag_of_Cameroon.png"), GL_CLAMP_TO_BORDER, GL_CLAMP_TO_BORDER);

}
void initialize_flag2(vec3 position_flag) {
	mesh_flag_support = mesh_primitive_cylinder(0.02f, { 0,0,0 }, { 0,0,3 });
	mesh_flag_support.push_back(mesh_primitive_cylinder(0.02f, { 0,0,3 }, { 0,1.0,3 }));
	flag_support = mesh_drawable(mesh_flag_support);


	mesh_flag = mesh_primitive_grid();
	flag2 = mesh_drawable(mesh_flag);
	flag2.transform.translate = position_flag;
	N_flag = sqrt(mesh_flag.position.size());
	for (int i = 0; i < N_flag; i++) {
		std::vector<vec3> temp2;
		velocities.push_back(temp2);
		for (int j = 0; j < N_flag; j++) {
			velocities[i].push_back({ 0,0,0 });
		}
	}

	L0 = 0.1f; // Rest length between A and B
	L0_cisaillement = L0 * sqrt(2);
	L0_flexion = 2 * L0;

	flag2.texture = opengl_texture_to_gpu(image_load_png("Assets/Village/Textures/Flag_of_CI.png"), GL_CLAMP_TO_BORDER, GL_CLAMP_TO_BORDER);

}

//rendering element to reflect or refract on corresponding frame buffers
void reflect_refract_to_fbos(WaterFrameBuffers* fbos, scene_environment& scene, GLFWwindow* window) {

	scene.clip_plan1 = vec4(0, 0, 0, 0);
	scene.clip_plan2 = vec4(0, 0, 0, 0);

	

	//storing current camera infos
	vec3 stored_camera_position = scene.camera.position();
	vec3 stored_camera_front = scene.camera.front();
	vec3 stored_camera_up = scene.camera.up();

	// rendering skybox, bridge and tree on reflection frame buffer
	fbos->bindReflectionFrameBuffer();

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);

	////Moving camera to miror the scene
	float dist = 2 * (scene.camera.position().z - water_plane.transform.translate.z);
	vec3 temp_pos = vec3(stored_camera_position[0], stored_camera_position[1], stored_camera_position[2] - dist);
	vec3 temp_front = vec3(stored_camera_front[0], stored_camera_front[1], -stored_camera_front[2]);
	vec3 temp_up = vec3(stored_camera_up[0], stored_camera_up[1], stored_camera_up[2]);

	scene.camera.look_at(temp_pos, temp_pos + temp_front, temp_up);
	////

	glDepthMask(GL_FALSE);
	draw_with_cubemap(skybox_drawable, scene);
	glDepthMask(GL_TRUE);
	create_full_ground_bridge(6, scene);

	//setting clipping planes to have selective rendering
	scene.clip_plan = vec4(0, 0, 1, -water_plane.transform.translate.z);
	
	tree1.transform.translate = ground.transform.translate;
	tree1.transform.translate.y += 5;
	tree1.transform.scale = 1.5;
	draw(tree1, scene);

	scene.clip_plan = vec4(0, 0, 0, 0);
	
	//considering fps camera movement
	if (camera_fps_active) {
		float const dt = timer.update();
		float moving_speed = 2.0f;
		//	scene.camera.position_camera += user.speed * 0.1f * dt * scene.camera.front();
		if (user.keyboard_state.up)
			cameraPosition += moving_speed * cameraFront * dt;
		//scene.camera.manipulator_translate_in_plane();
		if (user.keyboard_state.down)
			cameraPosition -= moving_speed * cameraFront * dt;
		update_camera();

		scene_fps.camera.look_at(cameraPosition, cameraPosition + cameraFront, cameraUp);

		if (user.keyboard_state.right)
			scene_fps.camera.manipulator_translate_in_plane(vec2(-dt * moving_speed, 0));
		if (user.keyboard_state.left)
			scene_fps.camera.manipulator_translate_in_plane(vec2(dt * moving_speed, 0));
		//scene.camera.manipulator_translate_in_plane();

		cameraPosition = scene_fps.camera.position();
		update_camera();
		scene_fps.camera.look_at(cameraPosition, cameraPosition + cameraFront, cameraUp);
		gun.transform.translate = cameraPosition + 0.3f * scene_fps.camera.right() - 0.15f * scene_fps.camera.up() + 0.7f * cameraFront;

		gun.transform.rotate = scene_fps.camera.orientation() * rotation({ 0,1,0 }, -pi / 2);
	}
	
	fbos->unbindCurrentFrameBuffer(window);

	
	// rendering ground on refraction frame buffer
	fbos->bindRefractionFrameBuffer();

	//setting camera to original coordinates
	scene.camera.look_at(stored_camera_position, stored_camera_position + stored_camera_front, stored_camera_up);

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);

	scene.clip_plan = vec4(0, 0, -1, water_plane.transform.translate.z);

	glDepthMask(GL_FALSE);
	draw_with_cubemap(skybox_drawable, scene);
	
	glDepthMask(GL_TRUE);

	
	ground.transform.scale = 30;
	ground.transform.translate = { 20, 0, -2.0f };
	draw(ground, scene);
	float dtf1 = timer.scale * 0.01;
	float mu1 = 0.002f;
	std::vector<std::vector<vcl::vec3>> temp;
	temp = animate_fish_hierarchy(fish_hierarchy, fish_positions, fish_velocities, fish_forces, m_fish, mu1, dtf1, r_fish, V_fish, 2 * pi);
	
	fish_positions = temp[0];
	fish_velocities = temp[1];
	fish_forces = temp[2];

	draw(tree1, scene);

	scene.clip_plan = vec4(0, 0, 0, 0);
	
	fbos->unbindCurrentFrameBuffer(window);
	

}

//putting fbos textures on water plane
void Send_Watertexture(WaterFrameBuffers* fbos, scene_environment& scene) {


	water_plane.texture = fbos->getReflectionTexture();
	water_plane.texture_2 = fbos->getRefractionTexture();
	water_plane.texture_3 = opengl_texture_to_gpu(image_load_png("Assets/water_assets/waterDUDV.png"), GL_REPEAT /**GL_TEXTURE_WRAP_S*/, GL_REPEAT /**GL_TEXTURE_WRAP_T*/);
	water_plane.texture_4 = opengl_texture_to_gpu(image_load_png("Assets/water_assets/water_normal.png"), GL_REPEAT /**GL_TEXTURE_WRAP_S*/, GL_REPEAT /**GL_TEXTURE_WRAP_T*/);
	water_plane.texture_5 = fbos->getRefractionDepthTexture();

	scene.cam_Pos = scene.camera.position();
	scene.t = timer1.update();
	Movefact += Wavespeed * scene.t;
	Movefact -= (int)Movefact;

}

mesh_drawable create_mouse_pointer() {
	mesh  pointer;
	pointer.push_back(mesh_primitive_triangle({ 0,0,0.001 }, { 0,0.001f, 0.002f }, { 0,-0.001f, 0.002f }));
	pointer.push_back(mesh_primitive_triangle({ 0,0,-0.001 }, { 0,0.001f, -0.002f }, { 0,-0.001f, -0.002f }));
	pointer.push_back(mesh_primitive_triangle({ 0,0.001,0 }, { 0,0.002f, 0.001f }, { 0,0.002f, -0.001f }));
	pointer.push_back(mesh_primitive_triangle({ 0,-0.001,0 }, { 0,-0.002f, 0.001f }, { 0,-0.002f, -0.001f }));
	return mesh_drawable(pointer);

}

int main(int, char* argv[])
{
	std::cout << "Run " << argv[0] << std::endl;

	int const width = 1280, height = 1024;
	GLFWwindow* window = create_window(width, height);
	window_size_callback(window, width, height);
	std::cout << opengl_info_display() << std::endl;;

	imgui_init(window);
	glfwSetCursorPosCallback(window, mouse_move_callback);
	glfwSetKeyCallback(window, keyboard_callback);
	glfwSetWindowSizeCallback(window, window_size_callback);

	WaterFrameBuffers* fbos = new WaterFrameBuffers(window);

	std::cout << "Initialize data ..." << std::endl;
	initialize_data();



	std::cout << "Start animation loop ..." << std::endl;
	user.fps_record.start();
	glEnable(GL_DEPTH_TEST);
	while (!glfwWindowShouldClose(window))
	{	
		glEnable(GL_CLIP_DISTANCE0);
		glEnable(GL_CLIP_DISTANCE1);
		glEnable(GL_CLIP_DISTANCE2);

		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		glClear(GL_DEPTH_BUFFER_BIT);

		if (camera_fps_active)
			scene_to_draw = &scene_fps;
		else
			scene_to_draw = &scene;
		
		(*scene_to_draw).light = (*scene_to_draw).camera.position();
		user.fps_record.update();

		
		imgui_create_frame();
		if (user.fps_record.event) {
			std::string const title = "VCL Display - " + str(user.fps_record.fps) + " fps";
			glfwSetWindowTitle(window, title.c_str());
		}
		

		ImGui::Begin("GUI", NULL, ImGuiWindowFlags_AlwaysAutoResize);
		user.cursor_on_gui = ImGui::IsAnyWindowFocused();

		if (user.gui.display_frame) draw(user.global_frame, *scene_to_draw);

		reflect_refract_to_fbos(fbos, *scene_to_draw, window);
		
		
		Send_Watertexture(fbos,*scene_to_draw);

		
		display_interface();
		display_frame();
		

		ImGui::End();
		imgui_render_frame(window);
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	fbos->cleanUp();
	imgui_cleanup();
	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}



void initialize_data()
{
	GLuint const shader_mesh = opengl_create_shader_program(read_text_file("shader/mesh.vert.glsl"), read_text_file("shader/mesh.frag.glsl"));
	GLuint const shader_uniform_color = opengl_create_shader_program(opengl_shader_preset("single_color_vertex"), opengl_shader_preset("single_color_fragment"));
	GLuint const shader = opengl_create_shader_program(read_text_file("shader/terrain_vert.glsl"), read_text_file("shader/terrain_frag.glsl"));
	GLuint const shader2 = opengl_create_shader_program(read_text_file("shader/terrain_main_vert.glsl"), read_text_file("shader/terrain_main_frag.glsl"));
	GLuint const shader_with_transparency = opengl_create_shader_program(read_text_file("shader/transparency_vert.glsl"), read_text_file("shader/transparency_frag.glsl"));
	GLuint const shader_skybox = opengl_create_shader_program(read_text_file("shader/skybox_vert.glsl"), read_text_file("shader/skybox_frag.glsl"));
	GLuint const shader_environment_map = opengl_create_shader_program(read_text_file("shader/environment_map_vert.glsl"), read_text_file("shader/environment_map_frag.glsl"));
	GLuint const shader_deform = opengl_create_shader_program(read_text_file("shader/shader_deform.vert.glsl"), read_text_file("shader/shader_deform.frag.glsl")); 

	GLuint const texture_white = opengl_texture_to_gpu(image_raw{ 1,1,image_color_type::rgba,{255,255,255,255} });

	mesh_drawable::default_shader = shader_mesh;
	mesh_drawable::default_texture = texture_white;
	curve_drawable::default_shader = shader_uniform_color;
	segments_drawable::default_shader = shader_uniform_color;


	user.global_frame = mesh_drawable(mesh_primitive_frame());
	user.gui.display_frame = false;
	scene.camera.distance_to_center = 2.5f;
	scene.camera.look_at({ 20,3,100 }, { 0,0,0 }, { 0,0,1 });

	scene_to_draw = &scene;

	parameters_central.frequency_gain = 2.47f;
	parameters_central.persistency = 0.55f;
	parameters_central.octave = 8;
	parameters_central.terrain_height = 0.36f;

	parameters_forest.frequency_gain = 2.53f;
	parameters_forest.persistency = 0.6f;
	parameters_forest.octave = 8;
	parameters_forest.terrain_height = 1.5f;



	// Read cubemap texture
	GLuint texture_cubemap = cubemap_texture();

	// Cube used to display the skybox
	skybox_drawable = mesh_drawable(mesh_primitive_cube({ 0,0,0 }, 50.0f), shader_skybox, texture_cubemap);
	



	image_raw const im_snow1 = image_load_png("Assets/snow/Textures/snow_10_diffuse.png");
	image_raw const im_snow2 = image_load_png("Assets/Terrain Textures Snow/TerrainTex_Snow2_07_d.png");
	image_raw const im_snow3 = image_load_png("Assets/snow/Textures/snow_02_diffuse.png");

	GLuint const texture_image_snow1_id = opengl_texture_to_gpu(im_snow1,
		GL_REPEAT /**GL_TEXTURE_WRAP_S*/,
		GL_REPEAT /**GL_TEXTURE_WRAP_T*/);
	GLuint const texture_image_snow2_id = opengl_texture_to_gpu(im_snow2,
		GL_REPEAT /**GL_TEXTURE_WRAP_S*/,
		GL_REPEAT /**GL_TEXTURE_WRAP_T*/);
	GLuint const texture_image_snow3_id = opengl_texture_to_gpu(im_snow3,
		GL_REPEAT /**GL_TEXTURE_WRAP_S*/,
		GL_REPEAT /**GL_TEXTURE_WRAP_T*/);


	image_raw const im_snow_path = image_load_png("Assets/Terrain Textures Snow/TerrainTex_Snow2_07_d.png");

	GLuint const texture_image_snow_path_id = opengl_texture_to_gpu(im_snow_path,
		GL_REPEAT /**GL_TEXTURE_WRAP_S*/,
		GL_REPEAT /**GL_TEXTURE_WRAP_T*/);

	// Create visual terrain surface
	std::vector<vcl::mesh> temp_terrains = create_main_terrain();
	std::vector<vcl::mesh> temp_terrains2 = create_mountain();
	mesh_forest = create_forest();
	mesh_main_terrain = temp_terrains[0];
	mesh_path_in_snow = temp_terrains[1];
	mesh_mountain_pos = temp_terrains2[0];
	mesh_mountain_neg = temp_terrains2[1];



	for (int i = 0; i < 9; i++) {
		positions_animals.push_back(evaluate_main_terrain(0.05f*i, 0.5f));
		positions_to_avoid.push_back(positions_animals[i]);
		velocities_animals.push_back({ 0,0,0 });
		destination_animals.push_back(evaluate_main_terrain(rand_interval(), rand_interval()));
		new_destination.push_back(false);
	}


	

	mouse_pointer = create_mouse_pointer();

	//creating water plane

	float scale = 1.0f;
	water_plane = mesh_drawable_multitexture3(mesh_primitive_grid({ -scale / 5,-scale,0 }, { scale / 5,-scale,0 }, { scale / 5,scale,0 }, { -scale / 5,scale,0 }, 10, 10));
	water_plane.shader = shader_deform;
	water_plane.transform.scale = 30;
	water_plane.transform.translate = { 20, 0, 1.2f };

	//creating water ground

	scale = 1.0f;
	ground = mesh_drawable(mesh_primitive_grid({ -scale,-scale,0 }, { scale,-scale,0 }, { scale,scale,0 }, { -scale,scale,0 }, 20, 20));
	ground.texture = opengl_texture_to_gpu(image_load_png("Assets/Village/Textures/Rock 1/Rock 1-diffuse.png"), GL_REPEAT, GL_REPEAT);



	//birds
	bird_hierarchy = create_bird_hierarchy();
	std::vector<std::vector<vec3>> temp_birds = initialize_birds(5, 10, -25, -3);
	std::vector<std::vector<vec3>> temp_birds1 = initialize_birds(4, 8,- 15.0, 5.0f);
	bird_positions = temp_birds[0];
	bird_positions.insert(bird_positions.end(), temp_birds1[0].begin(), temp_birds1[0].end());
	bird_velocities = temp_birds[1];
	bird_velocities.insert(bird_velocities.end(), temp_birds1[1].begin(), temp_birds1[1].end());
	bird_forces = temp_birds1[2];
	bird_forces.insert(bird_forces.end(), temp_birds1[2].begin(), temp_birds1[2].end());
	//end birds

	//fishs
	fish_hierarchy = create_fish_hierarchy();
	std::vector<std::vector<vec3>> temp_fishs = initialize_fishs(3, ground.transform.translate.z+0.1, water_plane.transform.translate.x-2 , water_plane.transform.translate.y-30 );
	std::vector<std::vector<vec3>> temp_fishs1 = initialize_fishs(4, ground.transform.translate.z + 0.1, water_plane.transform.translate.x+3, water_plane.transform.translate.y-25);
	fish_positions = temp_fishs[0];
	fish_positions.insert(fish_positions.end(), temp_fishs1[0].begin(), temp_fishs1[0].end());
	fish_velocities = temp_fishs[1];
	fish_velocities.insert(fish_velocities.end(), temp_fishs1[1].begin(), temp_fishs1[1].end());
	fish_forces = temp_fishs1[2];
	fish_forces.insert(fish_forces.end(), temp_fishs1[2].begin(), temp_fishs1[2].end());
	//end fishs

	//enemies
	enemy_hierarchy = create_bird_hierarchy();
	//end enemies

	//Uploading houses, and rock objects and textures
	mesh_houses.push_back(mesh_load_file_obj("Assets/Village/Models/Buildings/OBJbig_house_01.obj"));
	mesh_houses.push_back(mesh_load_file_obj("Assets/Village/Models/Buildings/OBJbuild_barracks_single_01.obj"));
	mesh_houses.push_back(mesh_load_file_obj("Assets/Village/Models/Buildings/OBJsmall_house_tall_roof_01.obj"));

	house_textures.push_back(opengl_texture_to_gpu(image_load_png("Assets/Village/Textures/Buildings/PNGbuild_building_01_a.png")));
	house_textures.push_back(opengl_texture_to_gpu(image_load_png("Assets/Village/Textures/Buildings/PNGbuild_building_01_a.png")));
	house_textures.push_back(opengl_texture_to_gpu(image_load_png("Assets/Village/Textures/Buildings/PNGbuild_building_01_a.png")));


	mesh_rock.push_back(mesh_load_file_obj("Assets/Rocks/Source/Models/OBJRock1_LP.obj"));
	mesh_rock.push_back(mesh_load_file_obj("Assets/Rocks/Source/Models/OBJRock1C.obj"));
	mesh_rock.push_back(mesh_load_file_obj("Assets/Rocks/Source/Models/OBJRock6C.obj"));


	rock_textures.push_back(opengl_texture_to_gpu(image_load_png("Assets/Rocks/Source/Textures/Rock1_snow.png")));
	rock_textures.push_back(opengl_texture_to_gpu(image_load_png("Assets/Rocks/Source/Textures/Rock1_snow.png")));
	rock_textures.push_back(opengl_texture_to_gpu(image_load_png("Assets/Rocks/Source/Textures/Rock6_snow.png")));

	


	//creating trees
	tree1 = mesh_drawable(create_tree());
	image_raw const im_trunk = image_load_png("Assets/Village/Textures/Vegetation/leaf_trunk.png");	
	GLuint const texture_image_trunk_id = opengl_texture_to_gpu(im_trunk,
		GL_REPEAT /**GL_TEXTURE_WRAP_S*/,
		GL_REPEAT /**GL_TEXTURE_WRAP_T*/);
	tree1.texture = texture_image_trunk_id;


	//generating the houses rocks and grass positions
	for (int i = 0; i < mesh_houses.size(); i++) {
		houses.push_back(mesh_drawable(mesh_houses[i]));
		houses[i].texture = house_textures[i];
		houses[i].shading.phong.specular = 0.0f;
		houses[i].shading.phong.diffuse = 0.2f;
	}


	houses_positions = generate_houses_positions(number_of_houses, 0.8f, 1.0f);
	rock_position = generate_houses_positions(number_of_rock, 0.8f, 0.5f);
	billboard_positions = generate_houses_positions(60, 0.3f, 0.5f);
	tree_position = generate_tree_positions(200, 1.0f, 1.0f);

	for (auto i : houses_positions) {
		house_indices.push_back((int)rand_interval(0.0f, houses.size()));
		houses_rand_rotation.push_back(rotation::axis_angle_to_matrix({ 0,0,1 }, rand_interval(0.0f, 2 * pi - 0.1f)));
		houses_rand_scale.push_back(rand_interval(0.03f, 0.05f));

	}


	for (int i = 0; i < mesh_rock.size(); i++) {
		rock.push_back(mesh_drawable(mesh_rock[i]));
		rock[i].texture = rock_textures[i];
		rock[i].shading.phong.specular = 0.0f;
		rock[i].shading.phong.diffuse = 0.2f;
	}

	for (auto i : rock_position) {
		rock_indices.push_back((int)rand_interval(0.0f, rock.size()));
		rock_rand_rotation.push_back(rotation::axis_angle_to_matrix({ 0,0,1 }, rand_interval(0.0f, 2 * pi - 0.1f)));
		rock_rand_scale.push_back(rand_interval(1.0f, 2.0f));

	}

	for (auto i : tree_position) {
		tree_rand_rotation.push_back(rotation::axis_angle_to_matrix({ 0,0,1 }, rand_interval(0.0f, 2 * pi - 0.1f)));
		tree_rand_scale.push_back(rand_interval(0.4f, 0.8f));
	}

	initial_house_rot = rotation::axis_angle_to_matrix(1 / sqrt(3) * vec3({ 1,1,1 }), 120.0f * pi / 180.0f);

	//construction of the bridge
	single_ground_plank = ground_plank();
	single_ground_transversal_plank = transversal_ground_plank();
	single_ground_longitudinal_plank = drawable_plank();


	forest_terrain = mesh_drawable_multitexture(mesh_forest);
	forest_terrain.shading.phong.specular = 0.0f; // non-specular terrain material
	forest_terrain.shading.color = { 1.0f,1.0f,1.0f };
	forest_terrain.texture = opengl_texture_to_gpu(image_load_png("Assets/Village/Textures/Grass 1/Grass 1-diffuse.png"), GL_REPEAT, GL_REPEAT);
	forest_terrain.texture_2 = opengl_texture_to_gpu(image_load_png("Assets/Village/Textures/Rock 1/Rock 1-diffuse.png"), GL_REPEAT, GL_REPEAT);
	forest_terrain.texture_3 = opengl_texture_to_gpu(image_load_png("Assets/Village/Textures/Grass 1/Grass 1-diffuse.png"), GL_REPEAT, GL_REPEAT);
	forest_terrain.shader = shader;


	main_terrain = mesh_drawable_multitexture(mesh_main_terrain);
	main_terrain.shading.phong.specular = 0.0f; // non-specular terrain material
	main_terrain.shading.color = { 1.0f,1.0f,1.0f };
	main_terrain.texture = texture_image_snow3_id;
	main_terrain.texture_2 = texture_image_snow1_id;
	main_terrain.texture_3 = texture_image_snow3_id;
	main_terrain.shader = shader;

	path_in_snow = mesh_drawable(mesh_path_in_snow);
	path_in_snow.shading.phong.specular = 0.0f; // non-specular terrain material
	path_in_snow.shading.color = { 1.0f,1.0f,1.0f };
	path_in_snow.texture = texture_image_snow_path_id;

	mountain_pos = mesh_drawable_multitexture(mesh_mountain_pos);
	mountain_pos.shading.phong.specular = 0.0f; // non-specular terrain material
	mountain_pos.shading.color = { 1.0f,1.0f,1.0f };
	mountain_pos.texture = texture_image_snow1_id;
	mountain_pos.texture_2 = texture_image_snow2_id;
	mountain_pos.texture_3 = texture_image_snow3_id;
	mountain_pos.shader = shader;


	mountain_neg = mesh_drawable_multitexture(mesh_mountain_neg);
	mountain_neg.shading.phong.specular = 0.0f; // non-specular terrain material
	mountain_neg.shading.color = { 1.0f,1.0f,1.0f };
	mountain_neg.texture = texture_image_snow1_id;
	mountain_neg.texture_2 = texture_image_snow2_id;
	mountain_neg.texture_3 = texture_image_snow3_id;
	mountain_neg.shader = shader;
	//main_terrain.shader = shader;



	mesh_billboard_grass = mesh_primitive_quadrangle({ 0, -0.5f,0 }, { 0, 0.5, 0 }, { 0.0, 0.5, 1.0 }, { 0, -0.5, 1.0f });
	billboard_grass = mesh_drawable(mesh_billboard_grass);
	billboard_grass.texture = opengl_texture_to_gpu(image_load_png("Assets/Village/Textures/grass.png"), GL_CLAMP_TO_BORDER, GL_CLAMP_TO_BORDER);
	billboard_grass.shading.phong.specular = 0.0f;
	billboard_grass.shader = shader_with_transparency;


	memory1 = vec3(mesh_billboard_grass.position[2].x, mesh_billboard_grass.position[2].y, mesh_billboard_grass.position[2].z);
	memory2 = vec3(mesh_billboard_grass.position[3].x, mesh_billboard_grass.position[3].y, mesh_billboard_grass.position[3].z);

	billboard_positions_grass = generate_grass(number_of_grass, 2.0f);
	for (int i = 0; i < billboard_positions_grass.size(); i++) {
		billboard_scale_grass.push_back(rand_interval(0.9f, 1.0f));
	}

	//set-up the two flags on both sides of the terrain
	initialize_flag({ 11.0f, -2, 4 });
	initialize_flag2({ 27.0f, -2, 4 });

	//Initialize snow
	float const r = 0.05f; // radius of the sphere
	snow_sphere = mesh_drawable(mesh_primitive_sphere(r));
	snow_sphere.shading.color = { 1.0f,1.0f,1.0f };
	
	//grass set-up
	mesh_billboard = mesh_primitive_quadrangle({ 0, -0.5f,0 }, { 0, 0.5, 0 }, { 0.0, 0.5, 1 }, { 0, -0.5, 1 });
	billboard = mesh_drawable(mesh_billboard);
	billboard.transform.scale = 0.5f;
	billboard.texture = opengl_texture_to_gpu(image_load_png("Assets/Village/Textures/Vegetation/grass_veg.png"), GL_CLAMP_TO_BORDER, GL_CLAMP_TO_BORDER);
	billboard.shader = shader_with_transparency;

	//used to animate the grass billboards. We store their positions to do an interpolation between their initial and their final position
	memory1 = vec3(mesh_billboard_grass.position[2].x, mesh_billboard_grass.position[2].y, mesh_billboard_grass.position[2].z);
	memory2 = vec3(mesh_billboard_grass.position[3].x, mesh_billboard_grass.position[3].y, mesh_billboard_grass.position[3].z);

	memory3 = vec3(mesh_billboard.position[2].x, mesh_billboard.position[2].y, mesh_billboard.position[2].z);
	memory4 = vec3(mesh_billboard.position[3].x, mesh_billboard.position[3].y, mesh_billboard.position[3].z);

	update_terrain(mesh_mountain_pos, mountain_pos, mesh_mountain_neg, mountain_neg, mesh_main_terrain, main_terrain, mesh_forest, forest_terrain, mesh_path_in_snow, path_in_snow, parameters, parameters_central, parameters_forest);
	update_houses_positions(houses_positions, mesh_main_terrain);
	update_houses_positions(vegetation_position, mesh_main_terrain);
	update_houses_positions(rock_position, mesh_main_terrain);
	update_houses_positions(billboard_positions, mesh_main_terrain);
	update_grass_positions(billboard_positions_grass, mesh_forest);
	update_tree_positions(tree_position, mesh_forest);
	update_houses_positions(destination_animals, mesh_main_terrain);
	positions_to_avoid.insert(positions_to_avoid.end(), houses_positions.begin(), houses_positions.end());
	positions_to_avoid.insert(positions_to_avoid.end(), rock_position.begin(), rock_position.end());

	init_camera_fps(mesh_main_terrain);
	
	//uploading the gun for the fps camera
	gun = mesh_drawable(mesh_load_file_obj("Assets/Gun/Models/OBJGun.obj"));
	gun.transform.scale = 15.0f;
	gun.transform.translate = cameraPosition + vec3(0.5f, -0.25f, -0.1f);
	gun.texture = opengl_texture_to_gpu(image_load_png("Assets//Gun/Textures/GunModle_Gun_AlbedoTransparency.png"));


	sphere = mesh_drawable(mesh_primitive_sphere(0.8f));
}

void update_camera() {
	if (cameraPosition[0] < 15.0f)
		update_camera_position(cameraPosition, mesh_main_terrain);
	if ( cameraPosition[0]>=15.0f && cameraPosition[0] < 25.0f)
		cameraPosition[2] = 2.3f;
	if (cameraPosition[0] >= 25.0f)
		update_camera_position(cameraPosition, mesh_forest);
}


void display_frame()
{
	//Switching between both cameras
	if (camera_fps_active)
		scene_to_draw = &scene_fps;
	else
		scene_to_draw = &scene;

	
	glDepthMask(GL_FALSE);
	draw_with_cubemap(skybox_drawable, *scene_to_draw);
	glDepthMask(GL_TRUE);

	float const dt = timer.update();
	float dtf = timer.scale * 0.01;
	float dt_billboard = timer_billboard.update();
	timer_enemy.update();
	t += dtf;
	

	//The timer to generate a new wind vector for the animation of the flags and the grass billboards
	bool const new_wind = timer_billboard.event;
	
	//Camera_fps
	if (camera_fps_active) {
		float moving_speed = 2.0f;
		//	scene.camera.position_camera += user.speed * 0.1f * dt * scene.camera.front();
		if (user.keyboard_state.up)
			cameraPosition += moving_speed * cameraFront * dt;
		//scene.camera.manipulator_translate_in_plane();
		if (user.keyboard_state.down)
			cameraPosition -= moving_speed * cameraFront * dt;
		update_camera();

		scene_fps.camera.look_at(cameraPosition, cameraPosition + cameraFront, cameraUp);
		
		if (user.keyboard_state.right)
			scene_fps.camera.manipulator_translate_in_plane(vec2(-dt * moving_speed, 0));
		if (user.keyboard_state.left)
			scene_fps.camera.manipulator_translate_in_plane(vec2(dt * moving_speed, 0));
		//scene.camera.manipulator_translate_in_plane();

		
		cameraPosition = scene_fps.camera.position();
		update_camera();
		scene_fps.camera.look_at(cameraPosition, cameraPosition + cameraFront, cameraUp);
		gun.transform.translate = cameraPosition + 0.3f * scene_fps.camera.right() - 0.15f * scene_fps.camera.up() + 0.7f * cameraFront;

		mouse_pointer.transform.translate = cameraPosition + 0.3f * cameraFront;
		mouse_pointer.transform.rotate = scene_fps.camera.orientation() * rotation({ 0,1,0 }, -pi / 2);

		gun.transform.rotate = scene_fps.camera.orientation() * rotation({ 0,1,0 }, -pi / 2);
		
		

		if (gun_animation) {
			
			gun.transform.translate = cameraPosition + 0.3f * scene_fps.camera.right() - 0.15f * scene_fps.camera.up() + 0.7f*fabs(cos(30*pi*(t-t_clicked))) * cameraFront;
			if (30 * pi * (t - t_clicked) >= pi)
				gun_animation = false;
		}
		draw(gun, *scene_to_draw);
		draw(mouse_pointer, *scene_to_draw);

	}
		
	
	draw(water_plane, *scene_to_draw); //drawing water plane 
	draw(main_terrain, mesh_main_terrain, *scene_to_draw);
	draw(path_in_snow, *scene_to_draw);
	draw(mountain_neg, mesh_mountain_neg, *scene_to_draw);
	draw(mountain_pos, mesh_mountain_pos, *scene_to_draw);
	draw(forest_terrain, mesh_forest, *scene_to_draw);
	create_full_ground_bridge(6, *scene_to_draw);
	
	//Populating our terrains with houses, rocks, tree...
	for (int i = 0; i < houses_positions.size(); i++) {


		houses[house_indices[i]].transform.translate = houses_positions[i];
		houses[house_indices[i]].transform.rotate = rotation(houses_rand_rotation[i]) * rotation(initial_house_rot);
		houses[house_indices[i]].transform.scale = houses_rand_scale[i];
		draw(houses[house_indices[i]], *scene_to_draw);
	}

	for (int i = 0; i < rock_position.size(); i++) {


		rock[rock_indices[i]].transform.translate = rock_position[i];
		rock[rock_indices[i]].transform.rotate = rotation(rock_rand_rotation[i]) * rotation(initial_house_rot);
		rock[rock_indices[i]].transform.scale = rock_rand_scale[i];
		draw(rock[rock_indices[i]], *scene_to_draw);
	}

	for (i = 0; i < tree_position.size(); i++) {

		tree1.transform.translate = tree_position[i];
		tree1.transform.translate.z -= 0.5;
		tree1.transform.rotate = rotation(tree_rand_rotation[i]);
		tree1.transform.scale = tree_rand_scale[i];
		draw(tree1, *scene_to_draw);
	}
	tree1.transform.translate = ground.transform.translate;
	tree1.transform.translate.y += 5;
	tree1.transform.scale = 1.5;
	
	draw(tree1, *scene_to_draw);
	
	

	//animation of the grass to which we apply a wind vector at a regular time interval
	if (new_wind) {
		last_timer = t;
		std::cout << last_timer << "\n";
		wind = vec3(rand_interval(-0.5f, 0.5f), rand_interval(-0.5f, 0.5f), 0.0f);
		wind_flag = vec3(rand_interval(-0.8f, 0.8f), 0.0f, 0.0f);
		friction = rand_interval(0.2f, 0.7f);
	}
	else
		wind_flag = vec3(0, 0, 0);
	
	mesh_billboard_grass.position[2] = memory1 + wind * std::exp(-friction * fabs(t - last_timer)) * std::sin(4 * (t - last_timer));
	mesh_billboard_grass.position[3] = memory2 - wind * std::exp(-friction * fabs(t - last_timer)) * std::sin(4 * (t - last_timer));


	mesh_billboard.position[2] = memory3 + wind * std::exp(-friction * fabs(t - last_timer)) * std::sin(4 * (t - last_timer));
	mesh_billboard.position[3] = memory4 - wind * std::exp(-friction * fabs(t - last_timer)) * std::sin(4 * (t - last_timer));

	billboard_grass.update_position(mesh_billboard_grass.position);
	billboard.update_position(mesh_billboard.position);

	for (int i = 0; i < billboard_positions_grass.size(); i++) {


		billboard_grass.transform.rotate = rotation();
		billboard_grass.transform.translate = billboard_positions_grass[i];
		billboard_grass.transform.scale = billboard_scale_grass[i];
		draw(billboard_grass, *scene_to_draw);
		billboard_grass.transform.rotate = rotation({ 0,0,1 }, pi / 2);
		billboard_grass.transform.translate = billboard_positions_grass[i];
		billboard_grass.transform.scale = billboard_scale_grass[i];
		draw(billboard_grass, *scene_to_draw);

	}

	//positioning the flag supports
	flag_support.transform.translate = {12, -2, 1};
	draw(flag_support, *scene_to_draw);
	flag_support.transform.translate = {28, -2, 1 };
	draw(flag_support, *scene_to_draw);


	float const m_flag = 0.01f;        // particle mass
	float const K = 10.0f;        // spring stiffness
	float mu = 0.002f;       // damping coefficient

	vec3 const g = { 0,0,-9.81f }; // gravity

	//animation of the birds
	std::vector<std::vector<vcl::vec3>> temp;
	
	temp = animate_bird_hierarchy(bird_hierarchy, bird_positions, bird_velocities, bird_forces, m_bird, mu, dtf, r_bird, V_bird, 2 * pi);
	bird_positions = temp[0];
	bird_velocities = temp[1];
	bird_forces = temp[2];
	//end_birds

	//animation of the enemies
	if (camera_fps_active) {
		std::vector<std::vector<vcl::vec3>> temp;
		if (timer_enemy.event) {
			temp = enemy_spawner(enemy_positions, enemy_velocities, enemy_forces, houses_positions, 1);
			enemy_positions = temp[0];
			enemy_velocities = temp[1];
			enemy_forces = temp[2];
		}
	
	
		temp = animate_enemy(enemy_hierarchy, enemy_positions, enemy_velocities, enemy_forces, m_bird, mu, dtf, r_bird, V_bird, 2 * pi);
		enemy_positions = temp[0];
		enemy_velocities = temp[1];
		enemy_forces = temp[2];
	}

	//end animation enemies
	
	//set up the forces in between the flag
	std::vector<std::vector<vec3>> forces;
	for (int i = 0; i < N_flag; i++) {
		std::vector<vec3> temp;
		forces.push_back(temp);
		for (int j = 0; j < N_flag; j++) {
			forces[i].push_back({ 0,0,0 });
		}

	}

	//Animation of the flag

	vec3 F_spring_tract1;
	vec3 F_spring_tract2;
	vec3 F_spring_tract3;
	vec3 F_spring_tract4;
	vec3 F_spring_cis1;
	vec3 F_spring_cis2;
	vec3 F_spring_cis3;
	vec3 F_spring_cis4;
	vec3 F_spring_flex1;
	vec3 F_spring_flex2;
	vec3 F_spring_flex3;
	vec3 F_spring_flex4;

		//Forces first

	for (unsigned int i = 0; i < N_flag - 1; i++) {
		for (unsigned int j = 0; j < N_flag; j++) {
			F_spring_tract1 = { 0,0,0 };
			F_spring_tract2 = { 0,0,0 };
			F_spring_tract3 = { 0,0,0 };
			F_spring_tract4 = { 0,0,0 };
			F_spring_cis1 = { 0,0,0 };
			F_spring_cis2 = { 0,0,0 };
			F_spring_cis3 = { 0,0,0 };
			F_spring_cis4 = { 0,0,0 };
			F_spring_flex1 = { 0,0,0 };
			F_spring_flex2 = { 0,0,0 };
			F_spring_flex3 = { 0,0,0 };
			F_spring_flex4 = { 0,0,0 };
			vec3 F_mass;
			vec3 F_friction;
			//ressorts de traction et cisaillement




			if (j != N_flag - 1) {
				F_spring_tract1 = spring_force(mesh_flag.position[j + i * N_flag], mesh_flag.position[i * N_flag + j + 1], L0, K);
				if (i != N_flag - 1) {
					F_spring_cis1 = spring_force(mesh_flag.position[i * N_flag + j], mesh_flag.position[(i + 1) * N_flag + j + 1], L0_cisaillement, K);
				}
			}
			if (i != N_flag - 1) {
				F_spring_tract2 = spring_force(mesh_flag.position[i * N_flag + j], mesh_flag.position[(i + 1) * N_flag + j], L0, K);
				if (j != 0) {
					F_spring_cis2 = spring_force(mesh_flag.position[i * N_flag + j], mesh_flag.position[(i + 1) * N_flag + j - 1], L0_cisaillement, K);
				}
			}

			if (j != 0) {
				F_spring_tract3 = spring_force(mesh_flag.position[i * N_flag + j], mesh_flag.position[i * N_flag + j - 1], L0, K);
				if (i != 0) {
					F_spring_cis3 = spring_force(mesh_flag.position[i * N_flag + j], mesh_flag.position[(i - 1) * N_flag + j - 1], L0_cisaillement, K);
				}
			}

			if (i != 0) {
				F_spring_tract4 = spring_force(mesh_flag.position[i * N_flag + j], mesh_flag.position[(i - 1) * N_flag + j], L0, K);

				if (j != N_flag - 1) {
					F_spring_cis4 = spring_force(mesh_flag.position[i * N_flag + j], mesh_flag.position[(i - 1) * N_flag + j + 1], L0_cisaillement, K);
				}
			}

			//ressorts de flexion
			if (j < N_flag - 2) {
				F_spring_flex1 = spring_force(mesh_flag.position[i * N_flag + j], mesh_flag.position[i * N_flag + j + 2], L0_flexion, K);
			}
			if (i < N_flag - 2) {
				F_spring_flex2 = spring_force(mesh_flag.position[i * N_flag + j], mesh_flag.position[(i + 2) * N_flag + j], L0_flexion, K);
			}
			if (j >= 2) {
				F_spring_flex3 = spring_force(mesh_flag.position[i * N_flag + j], mesh_flag.position[i * N_flag + j - 2], L0_flexion, K);
			}
			if (i >= 2) {
				F_spring_flex4 = spring_force(mesh_flag.position[i * N_flag + j], mesh_flag.position[(i - 2) * N_flag + j], L0_flexion, K);
			}

			//poids
			F_mass = m_flag * g;
			F_friction = -mu * velocities[i][j];

			forces[i][j] = F_spring_tract1 + F_spring_tract2 + F_spring_tract3 + F_spring_tract4 +
				F_spring_cis1 + F_spring_cis2 + F_spring_cis3 + F_spring_cis4 +
				F_spring_flex1 + F_spring_flex2 + F_spring_flex3 + F_spring_flex4 +
				F_mass + F_friction + wind_flag;



		}

	}


	for (int i = 0; i < N_flag - 1; i++) {
		for (int j = 0; j < N_flag; j++) {
			velocities[i][j] = (1 - mu) * velocities[i][j] + dtf * forces[i][j] / m_flag;
			mesh_flag.position[i * N_flag + j] = mesh_flag.position[i * N_flag + j] + dtf * velocities[i][j];

		}

	}

	flag.update_position(mesh_flag.position);
	draw(flag, *scene_to_draw);

	flag2.update_position(mesh_flag.position);
	draw(flag2, *scene_to_draw);

	
	//generation of the snow particules at different point of the field
//	bool const new_particle = timer.event;
//	if (new_particle == true) {
//		vec3 const p00 = { 0,0,20 };
//		vec3 const p01 = { 10,7,20 };
//		vec3 const p02 = { -10,7,20 };
//		vec3 const p03 = { -10,-7,20 };
//		vec3 const p04 = { 10,-7,20 };
//		vec3 const p05 = { 5,14,20 };
//		vec3 const p06 = { -5,14,20 };
//		vec3 const p07 = { -5,-14, 20 };
//		vec3 const p08 = { 5,-14,20 };
//		std::vector<vec3> vecp0 = { p00, p01, p02, p03, p04, p05, p06, p07, p08 };
		// Initial random velocity (x,y) components are uniformly distributed along a circle.



//		for (int i = 0; i < 9; i++) {
//			const float theta = rand_interval(0, 2 * pi);
//			const vec3 v0 = vec3(8.0f * std::cos(theta), 8.0f * std::sin(theta), -5.0f * std::sin(theta));
//
//			snow_particles.push_back({ vecp0[i],v0 });
//		}

//	}

	// Evolve position of particles
	
//	mu = 0.0025f;
//	for (particle_structure& particle : snow_particles)
//	{
//		const float m = 0.001f; // particle mass
//
//		vec3& p = particle.p;
//		vec3& v = particle.v;
//
//		const vec3 F = m * g - mu * v;
	
		// Numerical integration
//		v = v + dt * F / m;
//		p = p + dt * v;

//	}


	// Remove particles that are too low
//	for (auto it = snow_particles.begin(); it != snow_particles.end(); ) {
//		if (it->p.z < -3)
//			it = snow_particles.erase(it);
//		if (it != snow_particles.end())
//			++it;
//	}

	// Display particles
//	for (particle_structure& particle : snow_particles)
//	{
//		snow_sphere.transform.translate = particle.p;
//		draw(snow_sphere, *scene_to_draw);
//	}

	// Enable use of alpha component as color blending for transparent elements
//  new color = previous color + (1-alpha) current color
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Disable depth buffer writing
	//  - Transparent elements cannot use depth buffer
	//  - They are supposed to be display from furest to nearest elements


	for (int i = 0; i < billboard_positions.size(); i++) {
		billboard.transform.rotate = rotation();
		billboard.transform.translate = billboard_positions[i];
		draw(billboard, *scene_to_draw);
		billboard.transform.rotate = rotation({ 0,0,1 }, pi / 2);
		billboard.transform.translate = billboard_positions[i];
		draw(billboard, *scene_to_draw);
	}


	
}


void display_interface()
{
	ImGui::Checkbox("Frame", &user.gui.display_frame);
	ImGui::Checkbox("Wireframe", &user.gui.display_wireframe);
	ImGui::Checkbox("FPS Camera", &camera_fps_active);
}


void window_size_callback(GLFWwindow*, int width, int height)
{
	glViewport(0, 0, width, height);
	float const aspect = width / static_cast<float>(height);
	if(!camera_fps_active)
		scene.projection = projection_perspective(50.0f * pi / 180.0f, aspect, 0.1f, 1000.0f);
	
	if(camera_fps_active)
		scene_fps.projection = projection_perspective(50.0f * pi / 180.0f, aspect, 0.1f, 1000.0f);
}


void mouse_move_callback(GLFWwindow* window, double xpos, double ypos)
{
	if (!camera_fps_active) {
		vec2 const  p1 = glfw_get_mouse_cursor(window, xpos, ypos);
		vec2 const& p0 = user.mouse_prev;
		glfw_state state = glfw_current_state(window);

		auto& camera = scene.camera;
		if (!user.cursor_on_gui) {
			if (state.mouse_click_left && !state.key_ctrl)
				scene.camera.manipulator_rotate_trackball(p0, p1);
			if (state.mouse_click_left && state.key_ctrl)
				camera.manipulator_translate_in_plane(p1 - p0);
			if (state.mouse_click_right)
				camera.manipulator_scale_distance_to_center((p1 - p0).y);
		}

		user.mouse_prev = p1;

	}
	
	//Movement of the fps camera using spherical polar coordinates
	if (camera_fps_active) {
		if (first_mouse) {
			lastX = xpos;
			lastY = ypos;
			first_mouse = false;
		}

		float xoffset = (lastX - xpos);
		float yoffset = (lastY - ypos);
		lastX = xpos;
		lastY = ypos;

		//sensitivity of the camera
		float sensitivity = 0.2f;
		xoffset *= sensitivity;
		yoffset *= sensitivity;

		yaw += xoffset;
		pitch += yoffset;

		if (pitch > 89.0f)
			pitch = 89.0f;
		if (pitch < -89.0f)
			pitch = -89.0f;


		vec3 front;
		front.x = std::cos(yaw * pi / 180) * std::cos(pitch * pi / 180);
		front.z = std::sin(pitch * pi / 180);
		front.y = std::sin(yaw * pi / 180) * std::cos(pitch * pi / 180);


		cameraFront = front;
		std::cout << front;
		scene_fps.camera.look_at(cameraPosition, cameraPosition + cameraFront, cameraUp);

		//the gun follows the camera and is placed with x and z offsets.
		gun.transform.translate = cameraPosition + 0.3f * scene_fps.camera.right() - 0.15f * scene_fps.camera.up() + 0.7f * cameraFront;
		
		mouse_pointer.transform.translate = cameraPosition + 0.3f * cameraFront;

		glfw_state state = glfw_current_state(window);
		if (state.mouse_click_left) {
			user.mouse_click = true;
			gun_animation = true;
			t_clicked = t;
		}
			
		else
			user.mouse_click = false;
		


	}
}

void keyboard_callback(GLFWwindow*, int key, int, int action, int)
{
	if (key == GLFW_KEY_UP) {
		if (action == GLFW_PRESS) user.keyboard_state.up = true;
		if (action == GLFW_RELEASE) user.keyboard_state.up = false;
	}

	if (key == GLFW_KEY_DOWN) {
		if (action == GLFW_PRESS) user.keyboard_state.down = true;
		if (action == GLFW_RELEASE) user.keyboard_state.down = false;
	}

	if (key == GLFW_KEY_LEFT) {
		if (action == GLFW_PRESS) user.keyboard_state.left = true;
		if (action == GLFW_RELEASE) user.keyboard_state.left = false;
	}

	if (key == GLFW_KEY_RIGHT) {
		if (action == GLFW_PRESS) user.keyboard_state.right = true;
		if (action == GLFW_RELEASE) user.keyboard_state.right = false;
	}
}


void opengl_uniform(GLuint shader, scene_environment const& current_scene)
{
	opengl_uniform(shader, "projection", current_scene.projection);
	opengl_uniform(shader, "view", current_scene.camera.matrix_view());
	opengl_uniform(shader, "light", current_scene.light, false);
	opengl_uniform(shader, "time", current_scene.t, false);
	opengl_uniform(shader, "plane", current_scene.clip_plan, false);
	opengl_uniform(shader, "plane1", current_scene.clip_plan1, false);
	opengl_uniform(shader, "plane2", current_scene.clip_plan2, false);
	opengl_uniform(shader, "Movefactor", Movefact, false);
	opengl_uniform(shader, "camera_Pos", current_scene.cam_Pos, false);
}



