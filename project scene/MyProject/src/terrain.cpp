#include "terrain.hpp"
#include "mesh_drawable_multitexture.h"

using namespace vcl;

float prev_noise = 0.0f;
float y_max = 20.0f, x_max = 30.0f;
float y_max_village = 20.0f, x_max_village = 30.0f;
float y_max_for = 60.0f, x_max_for = 60.0f;

float path_width = 1.0f;


mesh initialize_second_terrain() {
	int const terrain_sample = 180;
	mesh terrain = mesh_primitive_grid({ -5,-10,0 }, { -25,-10,0 }, { -25,10,0 }, { -5,10,0 }, terrain_sample, terrain_sample);
	return terrain;
}

vec3 evaluate_mountain_pos(float u, float v) {

    float const x = x_max * (u - 0.5f);
    float const y = y_max * (v - 0.5f) ;


    std::vector<vec2> const ui = { {0.3f,0.4f}, {0.7f,0.6f}};
    std::vector<float> const hi = { 8.0f, 9.0f};
    std::vector<float> const sigmai = { 0.2f, 0.25f};


    int i = 0;
    float d = 0;
    float z_temp = 0;
    for (auto pi : ui) {
        d = norm(vec2(u, v) - pi) / sigmai[i];
        z_temp = z_temp + hi[i] * std::exp(-d * d);
        i++;
    }

    

    float const z = z_temp;

    return { x,y + y_max/2 + y_max_village/2-3.0f,z-0.1f};
}

vec3 evaluate_mountain_neg(float u, float v) {

    float const x = x_max * (u - 0.5f);
    float const y = y_max* (v - 0.5f) ;


    std::vector<vec2> const ui = {{0.2f,0.35f }, {0.55f,0.45f }, {0.7f, 0.4f} };
    std::vector<float> const hi = { 7.0f, 6.0f, 9.0f };
    std::vector<float> const sigmai = {0.25f, 0.2f, 0.2f };


    int i = 0;
    float d = 0;
    float z_temp = 0;

    
    for (auto pi : ui) {
        d = norm(vec2(u, v) - pi) / sigmai[i];
        z_temp = z_temp + hi[i] * std::exp(-d * d);
        i++;
    }

    float const z = z_temp;

    return { x,y - y_max / 2 - y_max_village / 2 +3.0f,z -0.4f};
}

vec3 evaluate_main_terrain(float u, float v)
{
    float const x = x_max_village * (u - 0.5f);
    float const y = y_max_village * (v - 0.5f);


    std::vector<vec2> const ui = { {0.1f,0.1f}, {0.1f,0.85f}, {0.6, 0.1f}, {0.6, 0.85f}};
    std::vector<float> const hi = { 2.0f, 3.0f, 1.5f, 2.5f};
    std::vector<float> const sigmai = { 0.25f, 0.35f, 0.25, 0.25};


    int i = 0;
    float d = 0;
    float z_temp = 0;

    for (auto pi : ui) {
       d = norm(vec2(u, v) - pi) / sigmai[i];
       z_temp = z_temp + hi[i] * std::exp(-d * d);
      i++;
     }    

    float const z = 1.2f + z_temp;

    return { x, y, z };
}

vec3 evaluate_forest_terrain(float u, float v)
{
    float const x = x_max_for * (u - 0.5f);
    float const y = y_max_for * (v - 0.5f);


    std::vector<vec2> const ui = { {0.8f,0.6f}, {0.8f,0.3f} , {0.8f, 0.2f}, {0.5f, 0.3f}, {0.2f, 0.8f}, {0.1f, 0.5f }, {0.1f, 0.8f} , {0.1f, 0.0f}, {0.1f, 0.88f}, {0.1f, 0.3f} };
    std::vector<float> const hi = { 4.0f, 5.0f, -1.0f, -3.0f, 5.0f , -1.0f, -2.0f, -1.0f, -1.5f, 0.5f};
    std::vector<float> const sigmai = { 0.2f, 0.25f, 0.3f, 0.4f, 0.35f, 0.35f , 0.35f, 0.35f, 0.35f, 0.35f};


    int i = 0;
    float d = 0;
    float z_temp = 0;
    for (auto pi : ui) {
        d = norm(vec2(u, v) - pi) / sigmai[i];
        z_temp = z_temp + hi[i] * std::exp(-d * d);
        i++;
    }



    float const z = z_temp;

    return { x + x_max/2 + x_max_for/2 + 10.0f, y, z };
}



vec3 evaluate_path_in_snow(float u, float v)
{
    float const x = x_max * (u - 0.5f);
    float const y = path_width * (v - 0.5f);

    //  vec2 const u0 = {0.5f, 0.5f};
     // float const h0 = 2.0f;
     // float const sigma0 = 0.15f;

    float const sigma = 0.0000025f;
    //float z_temp = (-10*std::exp(-0.0025f*y*y)  -std::exp(-0.0002f*y*y));
   // y = std::exp(-0.0025f * x * x);

    std::vector<vec2> const ui = { {0.3f,0.3f}, {0.5f,0.5f}, {0.2f, 0.7f}, {-0.3f, -1.0f },  {0.8f, 0.7f } };
    std::vector<float> const hi = { 9.0f, -5.0f , -5.0f , -5.0f, 5.0f };
    std::vector<float> const sigmai = { 0.3f, 0.5f, 0.5f, 0.5f, 0.1f };

    int i = 0;
    float d = 0;
    float z_temp = 2.0f;

//    for (auto pi : ui) {
 //       d = norm(vec2(u, v) - pi) / sigmai[i];
  //      z_temp = z_temp + hi[i] * std::exp(-d * d);
   //     i++;
    //}


    float const z = z_temp;


    return { x, y, z };
}



std::vector<vcl::mesh> create_mountain() {
    const unsigned int N = 50;

    mesh terrain_pos; // temporary terrain storage (CPU only)
    terrain_pos.position.resize(N * N);
    terrain_pos.uv.resize(N * N);

    mesh terrain_neg; // temporary terrain storage (CPU only)
    terrain_neg.position.resize(N * N);
    terrain_neg.uv.resize(N * N);
    

    // Fill terrain geometry
    for (unsigned int ku = 0; ku < N; ++ku)
    {
        for (unsigned int kv = 0; kv < N; ++kv)
        {
            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku / (N - 1.0f);
            const float v = kv / (N - 1.0f);

            // Compute the local surface function
            vec3 const p_pos = evaluate_mountain_pos(u, v);
            vec3 const p_neg = evaluate_mountain_neg(u, v);

            // Store vertex coordinates
            terrain_pos.position[kv + N * ku] = p_pos;
            terrain_pos.uv[kv + N * ku] = { 25 * u, 25 * v };

            terrain_neg.position[kv + N * ku] = p_neg;
            terrain_neg.uv[kv + N * ku] = { 25*u, 25*v };

        }
    }

    // Generate triangle organization
    //  Parametric surface with uniform grid sampling: generate 2 triangles for each grid cell
    for (size_t ku = 0; ku < N - 1; ++ku)
    {
        for (size_t kv = 0; kv < N - 1; ++kv)
        {
            const unsigned int idx = kv + N * ku; // current vertex offset

            const uint3 triangle_1 = { idx, idx + 1 + N, idx + 1 };
            const uint3 triangle_2 = { idx, idx + N, idx + 1 + N };

            terrain_pos.connectivity.push_back(triangle_1);
            terrain_pos.connectivity.push_back(triangle_2);
            terrain_neg.connectivity.push_back(triangle_1);
            terrain_neg.connectivity.push_back(triangle_2);
        }
    }


    terrain_pos.fill_empty_field(); // need to call this function to fill the other buffer with default values (normal, color, etc)
    terrain_neg.fill_empty_field();
    return { terrain_pos, terrain_neg };
}




std::vector<vcl::mesh> create_main_terrain() {
    const unsigned int N = 50;

    mesh terrain; // temporary terrain storage (CPU only)
    terrain.position.resize(N * N);
    terrain.uv.resize(N * N);
    mesh path_in_snow;
    unsigned int nu=0, nv=0;

    // Fill terrain geometry
    for (unsigned int ku = 0; ku < N; ++ku)
    {
        nu = nu + (unsigned int)1;
        nv = 0;
        for (unsigned int kv = 0; kv < N; ++kv)
        {
            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku / (N - 1.0f);
            const float v = kv / (N - 1.0f);

            // Compute the local surface function
            vec3 const p = evaluate_main_terrain(u, v);

            // Store vertex coordinates
            terrain.position[kv + N * ku] = p;
            terrain.uv[kv + N * ku] = { 15 * u, 15 * v };

            if (std::fabs(p[1]) <= path_width / 2) {
                path_in_snow.position.push_back(p);
                nv = nv + (unsigned int)1;
            }
        }
    }

    path_in_snow.uv.resize(nu * nv);
    for (unsigned int ku=0; ku < nu; ++ku) {
        for (unsigned int kv = 0; kv < nv; ++kv) {
            const float u = ku / (nu - 1.0f);
            const float v = kv / (nv - 1.0f);

            path_in_snow.uv[kv + nv * ku] = { 40 * u, 2 * v };
        }  
    }


    // Generate triangle organization
    //  Parametric surface with uniform grid sampling: generate 2 triangles for each grid cell
    for (size_t ku = 0; ku < N - 1; ++ku)
    {
        for (size_t kv = 0; kv < N - 1; ++kv)
        {
            const unsigned int idx = kv + N * ku; // current vertex offset

            const uint3 triangle_1 = { idx, idx + 1 + N, idx + 1 };
            const uint3 triangle_2 = { idx, idx + N, idx + 1 + N };

            terrain.connectivity.push_back(triangle_1);
            terrain.connectivity.push_back(triangle_2);
        }
    }

    for (size_t ku = 0; ku < nu - 1; ++ku)
    {
        for (size_t kv = 0; kv < nv - 1; ++kv)
        {
            const unsigned int idx = kv + nv * ku; // current vertex offset

            const uint3 triangle_1 = { idx, idx + 1 + nv, idx + 1 };
            const uint3 triangle_2 = { idx, idx + nv, idx + 1 + nv };

            path_in_snow.connectivity.push_back(triangle_1);
            path_in_snow.connectivity.push_back(triangle_2);
        }
    }

    

    terrain.fill_empty_field(); // need to call this function to fill the other buffer with default values (normal, color, etc)
    path_in_snow.fill_empty_field();
    return { terrain, path_in_snow };
}


vcl::mesh create_forest() {
    const unsigned int N = 50;

    mesh terrain; // temporary terrain storage (CPU only)
    terrain.position.resize(N * N);
    terrain.uv.resize(N * N);

    // Fill terrain geometry
    for (unsigned int ku = 0; ku < N; ++ku)
    {
        for (unsigned int kv = 0; kv < N; ++kv)
        {
            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku / (N - 1.0f);
            const float v = kv / (N - 1.0f);

            // Compute the local surface function
            vec3 const p = evaluate_forest_terrain(u, v);
            

            // Store vertex coordinates
            terrain.position[kv + N * ku] = p;
            terrain.uv[kv + N * ku] = { 25 * u, 25 * v };
        }
    }

    // Generate triangle organization
    //  Parametric surface with uniform grid sampling: generate 2 triangles for each grid cell
    for (size_t ku = 0; ku < N - 1; ++ku)
    {
        for (size_t kv = 0; kv < N - 1; ++kv)
        {
            const unsigned int idx = kv + N * ku; // current vertex offset

            const uint3 triangle_1 = { idx, idx + 1 + N, idx + 1 };
            const uint3 triangle_2 = { idx, idx + N, idx + 1 + N };

            terrain.connectivity.push_back(triangle_1);
            terrain.connectivity.push_back(triangle_2);
        }
    }


    terrain.fill_empty_field(); // need to call this function to fill the other buffer with default values (normal, color, etc)
    return terrain;
}



void update_terrain(mesh& mountain_pos, mesh_drawable& mountain_pos_visual,
         mesh& mountain_neg, mesh_drawable& mountain_neg_visual,
    mesh& main_terrain, mesh_drawable& main_terrain_visual,
    mesh& forest_terrain, mesh_drawable& forest_terrain_visual,
    mesh& path, mesh_drawable& path_visual,
    perlin_noise_parameters const& parameters, perlin_noise_parameters const& parameters_central, perlin_noise_parameters const& parameters_forest)
{
    // Number of samples in each direction (assuming a square grid)
    int const N_mountain = std::sqrt(mountain_pos.position.size());
    int const N_village = std::sqrt(main_terrain.position.size());
    int const N_forest = std::sqrt(forest_terrain.position.size());
    int idk = 0;
    int idv = 0;
    
    // Recompute the new vertices
    for (int ku = 0; ku < N_mountain; ++ku) {
        int temp_ku = 0;
        for (int kv = 0; kv < N_mountain; ++kv) {

            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku / (N_mountain - 1.0f);
            const float v = kv / (N_mountain - 1.0f);

            int const idx = ku * N_mountain + kv;

            // Compute the Perlin noise
            float const noise_mountain = noise_perlin({ u, v }, parameters.octave, parameters.persistency, parameters.frequency_gain);
           
            // use the noise as height value
            mountain_pos.position[idx].z =  evaluate_mountain_pos(u, v)[2] + parameters.terrain_height * noise_mountain;
            mountain_neg.position[idx].z = evaluate_mountain_neg(u, v)[2] + parameters.terrain_height * noise_mountain;

            // use also the noise as color value
            //terrain.color[idx] = noise * vec3(1, 1, 1); // à corriger plus tard.

 
        }
    }

    for (int ku = 0; ku < N_village; ku++) {
        int temp_ku = 0;
        for (int kv = 0; kv < N_village; kv++) {

            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku / (N_village - 1.0f);
            const float v = kv / (N_village - 1.0f);

            int const idx = ku * N_village + kv;

            // Compute the Perlin noise
            float const noise_main = noise_perlin({ u, v }, parameters_central.octave, parameters_central.persistency, parameters_central.frequency_gain);
           
            // use the noise as height value
            main_terrain.position[idx].z = evaluate_main_terrain(u, v)[2] + parameters_central.terrain_height * noise_main ;


            // use also the noise as color value
            //terrain.color[idx] = noise * vec3(1, 1, 1); // à corriger plus tard.

            if (std::fabs(main_terrain.position[idx].y) <= path_width / 2) {
                path.position[idk].z = main_terrain.position[idx].z + 0.001;
                path.color[idk] = 0.7 * noise_main * vec3(1, 1, 1);
                ++idk;
            }

        }
    }

    for (int ku = 0; ku < N_forest; ku++) {
        int temp_ku = 0;
        for (int kv = 0; kv < N_forest; kv++) {

            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku / (N_forest - 1.0f);
            const float v = kv / (N_forest - 1.0f);

            int const idx = ku * N_forest + kv;

            // Compute the Perlin noise
            float const noise_main = noise_perlin({ u, v }, parameters_forest.octave, parameters_forest.persistency, parameters_forest.frequency_gain);

            // use the noise as height value
            forest_terrain.position[idx].z = evaluate_forest_terrain(u, v)[2] + parameters_forest.terrain_height * noise_main;


            // use also the noise as color value
            //terrain.color[idx] = noise * vec3(1, 1, 1); // à corriger plus tard.

        }
    }

    // Update the normal of the mesh structure
    mountain_neg.compute_normal();
    mountain_pos.compute_normal();
    main_terrain.compute_normal();
    path.compute_normal();

    // Update step: Allows to update a mesh_drawable without creating a new one
    mountain_neg_visual.update_position(mountain_neg.position);
    mountain_neg_visual.update_normal(mountain_neg.normal);
    mountain_neg_visual.update_color(mountain_neg.color);

    mountain_pos_visual.update_position(mountain_pos.position);
    mountain_pos_visual.update_normal(mountain_pos.normal);
    mountain_pos_visual.update_color(mountain_pos.color);

    main_terrain_visual.update_position(main_terrain.position);
    main_terrain_visual.update_normal(main_terrain.normal);
    main_terrain_visual.update_color(main_terrain.color);
    
    forest_terrain_visual.update_position(forest_terrain.position);
    forest_terrain_visual.update_normal(forest_terrain.normal);
    forest_terrain_visual.update_color(forest_terrain.color);

    path_visual.update_position(path.position);
    path_visual.update_normal(path.normal);
    path_visual.update_color(path.color);

}

void update_houses_positions(std::vector<vec3> &houses_positions, mesh& terrain) {
    int const N = std::sqrt(terrain.position.size());
    float x, y, u, v;
    int ku, kv, idx;

    for (int i = 0; i < houses_positions.size(); i++) {
        x = houses_positions[i][0];
        y = houses_positions[i][1];

        u = (x)/x_max_village + 0.5f;
        v = y/y_max_village + 0.5f;
        ku = (int)((N - 1) * u);
        kv = (int)((N - 1) * v);
        idx = ku * N + kv;
        houses_positions[i][2] = terrain.position[idx].z;
    }
}

void update_camera_position(vec3& camera_position, mesh& terrain) {
    int const N = std::sqrt(terrain.position.size());
    float x, y, u, v, z_offset;
    int ku, kv, idx;
        
        x = camera_position[0];
        y = camera_position[1];
        if (x <= x_max_village/2){
            u = (x) / x_max_village + 0.5f;
            v = y / y_max_village + 0.5f;
            z_offset = 0.3f;
        }

        else{
            u = (x - x_max / 2 - x_max_for / 2 - 10.0f) / x_max_for + 0.5f;
            v = y / y_max_for + 0.5f;
            z_offset = 0.6f;
        }
        
        ku = (int)((N - 1) * u);
        kv = (int)((N - 1) * v);
        idx = ku * N + kv;
        camera_position[2] = terrain.position[idx].z + z_offset;

}


void update_rock(mesh& rock, mesh_drawable& rock_visual, perlin_noise_parameters const& parameters)
{
    // Number of samples in each direction (assuming a square grid)
    int const N = std::sqrt(rock.position.size());

    // Recompute the new vertices
    for (int ku = 0; ku < N; ++ku) {
        for (int kv = 0; kv < N; ++kv) {

            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku / (N - 1.0f);
            const float v = kv / (N - 1.0f);

            int const idx = ku * N + kv;

            // Compute the Perlin noise
            float const noise = noise_perlin({ u, v }, parameters.octave, parameters.persistency, parameters.frequency_gain);

            // use the noise as height value
            rock.position[idx].z = rock.position[idx].z + (parameters.terrain_height * noise-prev_noise);
            rock.position[idx].y = rock.position[idx].y + (parameters.terrain_height * noise-prev_noise);
            rock.position[idx].x = rock.position[idx].x + (parameters.terrain_height * noise- prev_noise);
            prev_noise = parameters.terrain_height * noise;
            // use also the noise as color value
           rock.color[idx] = 0.3f * vec3(0, 0.5f, 0) + 0.7f * noise * vec3(1, 1, 1); // à corriger plus tard.
        }
    }

    // Update the normal of the mesh structure
    rock.compute_normal();

    // Update step: Allows to update a mesh_drawable without creating a new one
    rock_visual.update_position(rock.position);
    rock_visual.update_normal(rock.normal);
    rock_visual.update_color(rock.color);

}

std::vector<vcl::vec3>  generate_houses_positions(int N, float safetyDistance, float safety_from_central_line) {
    std::vector<vec3> positions = {};
    vcl::vec3 eval;
    bool isIntersected;
    int j;
    float u, v;
    for (int i = 0; i < N; i++) {
        do {
            isIntersected = false;
            j = 0;
            do {
                u = rand_interval();
                v = rand_interval();
            } while (std::fabs(y_max_village * (v - 0.5f)) <= safety_from_central_line || 
                std::fabs(y_max_village * (v - 0.5f))>=y_max_village/2-2.0f || 
                std::fabs(x_max_village * (u - 0.5f)) >= x_max_village/2 - 2.0f
                );

            eval = evaluate_main_terrain(u, v);

            while (!isIntersected && j < i) {
                if (std::fabs(eval[0] - positions[j][0]) <= safetyDistance && std::fabs(eval[1] - positions[j][1]) <= safetyDistance) {
                    isIntersected = true;
                }
                j++;
            }

        } while (isIntersected);


        eval[2] += 0.0f;
        positions.push_back(eval);
    }

    return positions;
}

std::vector<vcl::vec3>  generate_grass(int N, float safety_from_central_line) {
    std::vector<vec3> positions = {};
    vcl::vec3 eval;
    bool isIntersected;
    int j;
    float u, v;
    float x_increment = x_max_for / N;
    float y_increment = y_max_for / N;
    float x = -x_max_for / 2 + x_increment, y;

    while (x < x_max_for / 2) {
        y = -y_max_for / 2 + y_increment;
        while (y < y_max_for / 2) {
            if (std::fabs(y) > safety_from_central_line) {
                u = (x) / x_max_for + 0.5f;
                v = y / y_max_for + 0.5f;
                positions.push_back(evaluate_forest_terrain(u, v));
            }
            y += y_increment;
        }
        
        x += x_increment;
    }

    return positions;
}

void update_grass_positions(std::vector<vec3>& grass_positions, mesh& terrain) {
    int const N = std::sqrt(terrain.position.size());
    float x, y, u, v;
    int ku, kv, idx;

    for (int i = 0; i < grass_positions.size(); i++) {
        x = grass_positions[i][0] -x_max / 2 - x_max_for / 2 - 10.0f;
        y = grass_positions[i][1];

        u = (x ) / x_max_for + 0.5f;
        v = y / y_max_for + 0.5f;
        ku = (int)((N - 1) * u);
        kv = (int)((N - 1) * v);
        idx = ku * N + kv;
        grass_positions[i][2] = terrain.position[idx].z;
    }
}

std::vector<vcl::vec3>  generate_tree_positions(int N, float safetyDistance, float safety_from_central_line) {
    std::vector<vec3> positions = {};
    vcl::vec3 eval;
    bool isIntersected;
    int j;
    float u, v;
    for (int i = 0; i < N; i++) {
        do {
            isIntersected = false;
            j = 0;
            do {
                u = rand_interval();
                v = rand_interval();
            } while (std::fabs(y_max_for * (v - 0.5f)) <= safety_from_central_line ||
                std::fabs(y_max_for * (v - 0.5f)) >= y_max_for / 2 - 1.0f ||
                std::fabs(x_max_for * (u - 0.5f)) >= x_max_for / 2 - 2.0f
                );

            eval = evaluate_forest_terrain(u, v);

            while (!isIntersected && j < i) {
                if (std::fabs(eval[0] - positions[j][0]) <= safetyDistance && std::fabs(eval[1] - positions[j][1]) <= safetyDistance) {
                    isIntersected = true;
                }
                j++;
            }

        } while (isIntersected);


        eval[2] += 0.0f;
        positions.push_back(eval);
    }

    return positions;
}

void update_tree_positions(std::vector<vec3>& tree_positions, mesh& terrain) {
    int const N = std::sqrt(terrain.position.size());
    float x, y, u, v;
    int ku, kv, idx;

    for (int i = 0; i < tree_positions.size(); i++) {
        x = tree_positions[i][0] - x_max / 2 - x_max_for / 2 - 10.0f;
        y = tree_positions[i][1];

        u = (x) / x_max_for + 0.5f;
        v = y / y_max_for + 0.5f;
        ku = (int)((N - 1) * u);
        kv = (int)((N - 1) * v);
        idx = ku * N + kv;
        tree_positions[i][2] = terrain.position[idx].z;
    }
}








