#pragma once

#include <iostream>
#include <string>
#include <vcl/vcl.hpp>
#include <GLFW/glfw3.h>
#include <gl/GLU.h>

static vcl::mesh loadOBJ(const char* file_name)
{
	//vertices
	std::vector<vcl::vec3> vertex_positions;
	std::vector<vcl::vec2> vertex_texcoords;
	std::vector<vcl::vec3> vertex_normals;

	//Face vectors
	std::vector<unsigned int> vertex_positions_indices;
	std::vector<std::vector<unsigned int>> vertex_texcoords_indices;
	std::vector<std::vector<unsigned int>> vertex_normals_indices;

	vcl::mesh mesh;

	std::stringstream ss;
	std::ifstream in_file(file_name);
	std::string line = "";
	std::string prefix = "";
	vcl::vec2 temp_vec2;
	vcl::vec3 temp_vec3;
	GLint temp_glint = 0;
	int count_indices = 0;

	if (!in_file.is_open()) {
		throw "ERROR:: OBJLOADER:: Could not open file!";
	}

	while (std::getline(in_file, line)) {
		//first the get prefix of the line
		ss.clear();
		ss.str(line);
		ss >> prefix;
		if (prefix == "#")
		{
			
		}
		else if (prefix == "o")
		{

		}
		else if (prefix == "s")
		{

		}
		else if (prefix == "use_mtl") 
		{

		}
		else if (prefix == "v")//vertices coordinates
		{
			ss >> temp_vec3.x >> temp_vec3.y >> temp_vec3.z;
			vertex_positions.push_back(temp_vec3);

		}
		else if (prefix == "vt")//vertices coordinates
		{
			ss >> temp_vec2.x >> temp_vec2.y;
			vertex_texcoords.push_back(temp_vec2);

		}
		else if (prefix == "vn")//vertices coordinates
		{
			ss >> temp_vec3.x >> temp_vec3.y >> temp_vec3.z;
			vertex_normals.push_back(temp_vec3);

		}
		else if (prefix == "f")//vertices coordinates
		{
			int counter = 0;
			int start = vertex_positions_indices.size();
			while (ss >> temp_glint) {
				if (counter == 0) 
				{
					vertex_positions_indices.push_back(temp_glint);
					++count_indices;

				}
				else if (counter == 1) 
				{
					vertex_texcoords_indices[vertex_positions_indices[vertex_positions_indices.size() - 1] - 1].push_back(temp_glint);
				}
				else if (counter == 2)
				{
					vertex_normals_indices[vertex_positions_indices[vertex_positions_indices.size()-1]-1].push_back(temp_glint);
				}

				if (ss.peek() == '/') 
				{
					++counter;
					ss.ignore(1, '/');
				}

				else if (ss.peek() == ' ') 
				{
					++counter;
					ss.ignore(1, ' ');
				}

				if (counter > 2)
				{
					counter = 0;
				}
			}

			if (count_indices == 3) {
				const vcl::uint3 triangle = {vertex_positions_indices[start]-1, vertex_positions_indices[start + 2]-1, vertex_positions_indices[start + 1] -1};
				mesh.connectivity.push_back(triangle);
				if (mesh.connectivity.size() ==  6408)
					std::cout << "problem here 1 " << triangle;
			}

			else if (count_indices > 3) 
			{
				for (int i = 0; i < count_indices-1; i += 3)
				{
					if (count_indices - i - 1 == 0) 
					{
					}

					else if(start + count_indices - i - 1 == 1)
					{
						const vcl::uint3 triangle = {vertex_positions_indices[start + i]-1, vertex_positions_indices[start + i + 2]-1, vertex_positions_indices[start + i + 1]-1 };
						mesh.connectivity.push_back(triangle);
						if (mesh.connectivity.size() == 6408)
							std::cout << "problem here 1 " << triangle;
						break;
					}
					
					else
					{
						const vcl::uint3 triangle_1 = {vertex_positions_indices[start + i]-1,
							vertex_positions_indices[start + i + 2] - 1, vertex_positions_indices[start + i + 1] - 1 };

						const vcl::uint3 triangle_2 = { vertex_positions_indices[start + i]-1,
						vertex_positions_indices[start + i + 3] - 1, vertex_positions_indices[start + i + 2] - 1 };

						mesh.connectivity.push_back(triangle_1);
						if (mesh.connectivity.size() == 6408)
							std::cout << "problem here 2 " << triangle_1;
						mesh.connectivity.push_back(triangle_2);
						if (mesh.connectivity.size() == 6408)
							std::cout << "problem here 3 " << triangle_2;

					}
					
				}

			}

			else {
				std::cout << "something wrong";
			}
			

			count_indices = 0;

		}
		else
		{
			
		}

		
		mesh.position.resize(vertex_positions.size());
		mesh.uv.resize(vertex_positions.size());
		mesh.normal.resize(vertex_positions.size());
		mesh.color.resize(vertex_positions.size());
		std::cout << "Nb vertices: " << mesh.position.size() << "\n";

		
	}
	
	for (size_t i = 0; i < mesh.position.size(); ++i) {
		vcl::vec3 normals_indices = { 0,0,0 };
		vcl::vec2 texcoords_indices = { 0,0};
		mesh.position[i] = vertex_positions[i];
		for (auto j : vertex_normals_indices[i]) {
			normals_indices += vertex_normals[j];	
		}
		for (auto j : vertex_texcoords_indices[i]) {
			texcoords_indices += vertex_texcoords[j]/vertex_texcoords_indices[i].size();
		}

		mesh.normal[i] = normals_indices / (std::sqrt(normals_indices[0] * normals_indices[0] + normals_indices[1] * normals_indices[1] +
			normals_indices[2] * normals_indices[2]));

		mesh.uv[i] = texcoords_indices;
		mesh.color[i] = vcl::vec3(1.0f, 1.0f, 1.0f);


	}

	std::cout << vertex_positions_indices[vertex_positions_indices.size() - 1];
//	for (size_t i = 0; i < mesh.position.size(); ++i) {
//		mesh.position[vertex_positions_indices[i]-1] = vertex_positions[vertex_positions_indices[i]-1];
//		mesh.uv[vertex_texcoords_indices[i]-1] = vertex_texcoords[vertex_texcoords_indices[i]-1];
//		mesh.normal[vertex_normals_indices[i]-1] = vertex_normals[vertex_normals_indices[i]-1];
//		mesh.color[vertex_positions_indices[i]-1] = vcl::vec3(1.0f,1.0f,1.0f);
//	}
	mesh.fill_empty_field();
	std::cout << "Nb vertices: " << mesh.position.size()<<"\n";
	std::cout << "OBJ file loaded!\n";

	return mesh;

}
