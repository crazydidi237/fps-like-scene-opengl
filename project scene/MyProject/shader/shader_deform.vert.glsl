#version 330 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec3 color;
layout (location = 3) in vec2 uv;

out struct fragment_data
{
    vec3 position;
    vec3 color;
	vec3 normal;
    vec2 uv;
	vec3 eye;
	vec4 gl_Position;
	vec3 ToCamVec;
} fragment;

const float tiling=1.0;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform vec3 camera_Pos;

uniform float time;

void main()
{
	vec4 worldPosition= model*vec4(position,1.0);
	fragment.ToCamVec=camera_Pos-worldPosition.xyz;
	vec3 p0 = position.xyz;
	
	vec3 p = vec3(p0.x, p0.y, 0.0 );
	// Compute normals after deformation
	
	
	
	fragment.position = vec3(model * vec4(p,1.0));
	fragment.color = color;
	fragment.uv = vec2(position.x/2.0 + 0.5, position.y/2.0 + 0.5)*tiling;
	//fragment.uv = vec2(position.x, position.y)*tiling;
	fragment.eye = vec3(inverse(view)*vec4(0,0,0,1.0));

	fragment.gl_Position =projection * view * model * vec4(p, 1.0);
	gl_Position =fragment.gl_Position;
	 
}
