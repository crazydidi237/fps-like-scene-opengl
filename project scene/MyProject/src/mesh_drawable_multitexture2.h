#pragma once

#include "vcl/vcl.hpp"

// Create a mesh_drawable_multitexture in reusing the mesh_drawable structure and adding an extra texture ID
struct mesh_drawable_multitexture2 : vcl::mesh_drawable {

	using mesh_drawable::mesh_drawable; // reuse the same constructor

	GLuint texture_2; // add a second texture

	GLuint texture_3;

	GLuint texture_4;
};


// Adapt the draw function for a mesh_drawable_multitexture
template <typename SCENE>
void draw(mesh_drawable_multitexture2 const& drawable, vcl::mesh const& m, SCENE const& scene)
{
	// Reuse the same draw function, but send an additional uniform for the second texture
	// ...

	assert_vcl(drawable.shader != 0, "Try to draw mesh_drawable without shader");
	assert_vcl(drawable.texture != 0, "Try to draw mesh_drawable without texture");
	glUseProgram(drawable.shader); opengl_check;
	opengl_uniform(drawable.shader, scene);
	opengl_uniform(drawable.shader, drawable.shading);
	opengl_uniform(drawable.shader, "model", drawable.transform.matrix());

	std::vector<float> altitudes;
	float max_height = -1000.0f, min_height = 1000.0f;
	float max_width = -1000.0f, min_width = 1000.0f;
	for (auto i : m.position) {
		if (i[2] > max_height) {
			max_height = i[2];
		}
		if (i[2] < min_height) {
			min_height = i[2];
		}
		if (i[1] > max_width) {
			max_width = i[1];
		}
		if (i[1] < min_width) {
			min_width = i[1];
		}
	}



	// Set the four textures
	glActiveTexture(GL_TEXTURE0); opengl_check;
	glBindTexture(GL_TEXTURE_2D, drawable.texture); opengl_check;
	glUniform1f(glGetUniformLocation(drawable.shader, "region1.max"), max_height); opengl_check;
	glUniform1f(glGetUniformLocation(drawable.shader, "region1.min"), min_height + 0.5f * (max_height - min_height + 0.5f)); opengl_check;
	glUniform1f(glGetUniformLocation(drawable.shader, "region1.ymax"),max_width); opengl_check;
	glUniform1f(glGetUniformLocation(drawable.shader, "region1.ymin"),0.05*max_width); opengl_check;
	vcl::opengl_uniform(drawable.shader, "region1ColorMap", 0);  opengl_check;

	glActiveTexture(GL_TEXTURE1); opengl_check; // the additional texture 
	glBindTexture(GL_TEXTURE_2D, drawable.texture_2); opengl_check;
	glUniform1f(glGetUniformLocation(drawable.shader, "region2.max"), min_height - 0.5f + 0.7f * (max_height - min_height + 0.5f)); opengl_check;
	glUniform1f(glGetUniformLocation(drawable.shader, "region2.min"), min_height - 0.5f + 0.05f * (max_height - min_height + 0.5f)); opengl_check;
	glUniform1f(glGetUniformLocation(drawable.shader, "region2.ymax"),0.1*max_width); opengl_check;
	glUniform1f(glGetUniformLocation(drawable.shader, "region2.ymin"),-0.1*max_width); opengl_check;
	vcl::opengl_uniform(drawable.shader, "region2ColorMap", 1);  opengl_check; // the second texture is called "image_texture_2" in the shader in this case

	glActiveTexture(GL_TEXTURE2); opengl_check; // the additional texture 
	glBindTexture(GL_TEXTURE_2D, drawable.texture_3); opengl_check;
	glUniform1f(glGetUniformLocation(drawable.shader, "region3.max"), min_height - 0.5f + 0.1f * (max_height - min_height + 0.5f)); opengl_check;
	glUniform1f(glGetUniformLocation(drawable.shader, "region3.min"), min_height - 0.5f); opengl_check;
	glUniform1f(glGetUniformLocation(drawable.shader, "region3.ymax"),-0.05* max_width); opengl_check;
	glUniform1f(glGetUniformLocation(drawable.shader, "region3.ymin"), min_width); opengl_check;
	vcl::opengl_uniform(drawable.shader, "region3ColorMap", 2);  opengl_check; // the second texture is called "image_texture_2" in the shader in this case


	// Standard call function
	assert_vcl(drawable.number_triangles > 0, "Try to draw mesh_drawable with 0 triangles"); opengl_check;
	glBindVertexArray(drawable.vao);   opengl_check;
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, drawable.vbo.at("index")); opengl_check;
	glDrawElements(GL_TRIANGLES, GLsizei(drawable.number_triangles * 3), GL_UNSIGNED_INT, nullptr); opengl_check;


	glBindVertexArray(0);

	// Clean the two textures binding to avoid any side effect after this draw
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, 0);

}