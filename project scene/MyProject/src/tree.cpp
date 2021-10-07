#include "tree.h"
#define _USE_MATH_DEFINES

#include <math.h>
using namespace vcl;



mesh create_tree_trunk_cylinder(float r, float h, vec3 Deviation, vec3 Pos_init)
{


    // Number of samples of the terrain is N x N
    const unsigned int N = 4;

    mesh cylinder; // temporary terrain storage (CPU only)
    cylinder.position.resize(N * N);
    cylinder.uv.resize(N * N);

    // Fill terrain geometry
    for (unsigned int ku = 0; ku < N; ++ku)
    {
        for (unsigned int kv = 0; kv < N; ++kv)
        {
            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku / (N - 1.0f);
            const float v = kv / (N - 1.0f);

            // Compute the local surface function
            vec3 const p = { r * std::cos(2 * pi * u) + Pos_init.x + v * Deviation.x, r * std::sin(2 * pi * u) + Pos_init.y + v * Deviation.y, h * (v)+Pos_init.z };
            vec2 const uv = { u,v };

            // Store vertex coordinates
            cylinder.position[kv + N * ku] = p;
            cylinder.uv[kv + N * ku] = { u / 2.,v / 2. };
        }
    }
    // Generate triangle organization
    for (size_t ku = 0; ku < N - 1; ++ku)
    {
        for (size_t kv = 0; kv < N - 1; ++kv)
        {
            const unsigned int idx = kv + N * ku;

            const uint3 triangle_1 = { idx, idx + 1 + N, idx + 1 };
            const uint3 triangle_2 = { idx, idx + N, idx + 1 + N };

            cylinder.connectivity.push_back(triangle_1);
            cylinder.connectivity.push_back(triangle_2);
        }
    }

    cylinder.fill_empty_field();
    return cylinder;
}

mesh create_leaf(vec3 Deviation, vec3 Pos_init)
{
    Deviation.z += 0.5;
    mesh leaf;
    leaf.position.resize(4);
    leaf.uv.resize(4);
    leaf.position[0] = Pos_init;
    leaf.position[1] = Pos_init + Deviation;
    vec3 normal;

    if (Deviation.y != 0) {
        normal.x = 1; normal.y = 0;
    }
    else {
        normal.x = 0; normal.y = 1;
    }

    leaf.position[2] = Pos_init + Deviation / 3 + normal / 4;
    leaf.position[3] = Pos_init + Deviation / 3 + normal / (-4);

    uint3 triangle_1 = { 0, 2, 3 };
    uint3 triangle_2 = { 3, 2, 1 };

    leaf.connectivity.push_back(triangle_1);
    leaf.connectivity.push_back(triangle_2);
    leaf.color.fill({ 0.f, 1.f, 0.f });
    leaf.uv[0] = { 0,0.5 };
    leaf.uv[1] = { 1,0.5 };
    leaf.uv[2] = { 0,1 };
    leaf.uv[3] = { 1,1 };


    leaf.fill_empty_field();
    return leaf;
}

void tree_rec(float r, float h, vec3 Deviation, vec3 Pos_init, mesh& Abr, int n) {
    if (n <= 1) {
        Abr.push_back(create_leaf(Deviation, Pos_init));
        return;
    }

    float proba = 0.8;
    if (rand_interval() <= proba) {
        float a = rand_interval(0.8, 1); float b = rand_interval(0.8, 1);  // a module le rayon des branches enfants, b la hauteur de ces branches  
        r = r * a;
        h = h * b;
        float c = rand_interval(0.0, 1); //c correspond à la pente crée par ces branche par rapport au parent
        Deviation.y = 1;
        Deviation *= h * c;
        vec3 ajust = { 0,r / a - r ,0 };
        vec3 Pos_ = Pos_init + ajust;
        Abr.push_back(create_tree_trunk_cylinder(r, h, Deviation, Pos_));
        vec3 Pos_new1 = Pos_ + Deviation;
        Pos_new1.z += h;
        tree_rec(r, h, Deviation, Pos_new1, Abr, n - 1);
    }

    if (rand_interval() <= proba) {
        float a = rand_interval(0.8, 1); float b = rand_interval(0.8, 1);   // a module le rayon des branches enfants, b la hauteur de ces branches  
        r = r * a;
        h = h * b;
        float c = rand_interval(0.3, 1);
        Deviation = { 0,0,0 };
        Deviation.x = -1;
        Deviation *= h * c;
        vec3 ajust = { -r / a + r,0,0 };
        vec3 Pos_ = Pos_init + ajust;
        Abr.push_back(create_tree_trunk_cylinder(r, h, Deviation, Pos_));
        vec3 Pos_new2 = Pos_ + Deviation;
        Pos_new2.z += h;
        tree_rec(r, h, Deviation, Pos_new2, Abr, n - 1);
    }

    if (rand_interval() <= proba) {
        float a = rand_interval(0.8, 1); float b = rand_interval(0.8, 1);   // a module le rayon des branches enfants, b la hauteur de ces branches  
        r = r * a;
        h = h * b;
        float c = rand_interval(0.3, 1);
        Deviation = { 0,0,0 };
        Deviation.y = -1;
        Deviation *= h * c;
        vec3 ajust = { 0,-r / a + r,0 };
        vec3 Pos_ = Pos_init + ajust;
        Abr.push_back(create_tree_trunk_cylinder(r, h, Deviation, Pos_));
        vec3 Pos_new3 = Pos_ + Deviation;
        Pos_new3.z += h;
        tree_rec(r, h, Deviation, Pos_new3, Abr, n - 1);
    }

    if (rand_interval() <= proba) {
        float a = rand_interval(0.8, 1); float b = rand_interval(0.8, 1);   // a module le rayon des branches enfants, b la hauteur de ces branches  
        r = r * a;
        h = h * b;
        float c = rand_interval(0.3, 1);
        Deviation = { 0,0,0 };
        Deviation.x = 1;
        Deviation *= h * c;
        vec3 ajust = { r / a - r,0,0 };
        vec3 Pos_ = Pos_init + ajust;
        Abr.push_back(create_tree_trunk_cylinder(r, h, Deviation, Pos_));
        vec3 Pos_new4 = Pos_ + Deviation;
        Pos_new4.z += h;
        tree_rec(r, h, Deviation, Pos_new4, Abr, n - 1);
    }

    if (rand_interval() <= proba) {
        float a = rand_interval(0.8, 1); float b = rand_interval(0.8, 1);   // a module le rayon des branches enfants, b la hauteur de ces branches  
        r = r * a;
        h = h * b;
        float c = rand_interval(0.3, 1);
        Deviation = { 0,0,0 };
        Deviation.x = 0.7;
        Deviation.y = 0.7;

        Deviation *= h * c;
        vec3 ajust = { 0,0,0 };
        vec3 Pos_ = Pos_init + ajust;
        Abr.push_back(create_tree_trunk_cylinder(r, h, Deviation, Pos_));
        vec3 Pos_new5 = Pos_ + Deviation;
        Pos_new5.z += h;
        tree_rec(r, h, Deviation, Pos_new5, Abr, n - 1);
    }

    if (rand_interval() <= proba) {
        float a = rand_interval(0.8, 1); float b = rand_interval(0.8, 1);   // a module le rayon des branches enfants, b la hauteur de ces branches  
        r = r * a;
        h = h * b;
        float c = rand_interval(0.3, 1);
        Deviation = { 0,0,0 };
        Deviation.x = -0.7;
        Deviation.y = -0.7;

        Deviation *= h * c;
        vec3 ajust = { 0,0,0 };
        vec3 Pos_ = Pos_init + ajust;
        Abr.push_back(create_tree_trunk_cylinder(r, h, Deviation, Pos_));
        vec3 Pos_new6 = Pos_ + Deviation;
        Pos_new6.z += h;
        tree_rec(r, h, Deviation, Pos_new6, Abr, n - 1);
    }
    if (rand_interval() <= proba) {
        float a = rand_interval(0.8, 1); float b = rand_interval(0.8, 1);  // a module le rayon des branches enfants, b la hauteur de ces branches  
        r = r * a;
        h = h * b;
        float c = rand_interval(0.3, 1);
        Deviation = { 0,0,0 };
        Deviation.x = 0.7;
        Deviation.y = -0.7;

        Deviation *= h * c;
        vec3 ajust = { 0,0,0 };
        vec3 Pos_ = Pos_init + ajust;
        Abr.push_back(create_tree_trunk_cylinder(r, h, Deviation, Pos_));
        vec3 Pos_new7 = Pos_ + Deviation;
        Pos_new7.z += h;
        tree_rec(r, h, Deviation, Pos_new7, Abr, n - 1);
    }
    if (rand_interval() <= proba) {
        float a = rand_interval(0.8, 1); float b = rand_interval(0.8, 1);   // a module le rayon des branches enfants, b la hauteur de ces branches  
        r = r * a;
        h = h * b;
        float c = rand_interval(0.3, 1);
        Deviation = { 0,0,0 };
        Deviation.x = -0.7;
        Deviation.y = 0.7;

        Deviation *= h * c;
        vec3 ajust = { 0,0,0 };
        vec3 Pos_ = Pos_init + ajust;
        Abr.push_back(create_tree_trunk_cylinder(r, h, Deviation, Pos_));
        vec3 Pos_new8 = Pos_ + Deviation;
        Pos_new8.z += h;
        tree_rec(r, h, Deviation, Pos_new8, Abr, n - 1);
    }

}


mesh create_tree()
{
    float h = 1; float r = h / 15.; vec3 Pos_init = { 0,0,h / 2.0 }; vec3 Deviation = { 0,0,0 }; int n = 5;
    mesh Abr = create_tree_trunk_cylinder(r, h, { 0,0,0 }, Pos_init);
    Pos_init.z += h;
    tree_rec(r, h, Deviation, Pos_init, Abr, n);
    return Abr;
}
