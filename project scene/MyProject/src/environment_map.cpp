#include "environment_map.hpp"

using namespace vcl;

GLuint cubemap_texture()
{
	// Load images
	image_raw const left = image_load_png("Assets/SkyboxSet1/CloudyLightRays/CloudyLightRaysLeft2048.png");
	image_raw const right = image_load_png("Assets/SkyboxSet1/CloudyLightRays/CloudyLightRaysRight2048.png");
	image_raw const top = image_load_png("Assets/SkyboxSet1/CloudyLightRays/CloudyLightRaysUp2048.png");
	image_raw const bottom = image_load_png("Assets/SkyboxSet1/CloudyLightRays/CloudyLightRaysDown2048.png");
	image_raw const front = image_load_png("Assets/SkyboxSet1/CloudyLightRays/CloudyLightRaysFront2048.png");
	image_raw const back = image_load_png("Assets/SkyboxSet1/CloudyLightRays/CloudyLightRaysBack2048.png");

	// Send images to GPU as cubemap
	GLuint cubemap;
	glGenTextures(1, &cubemap);
	glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap);

	glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_X, 0, GL_RGBA4, left.width, left.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, ptr(left.data));
	glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X, 0, GL_RGBA4, right.width, right.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, ptr(right.data));
	glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_Y, 0, GL_RGBA4, front.width, front.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, ptr(front.data));
	glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, 0, GL_RGBA4, back.width, back.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, ptr(back.data));
	glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_Z, 0, GL_RGBA4, top.width, top.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, ptr(top.data));
	glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, 0, GL_RGBA4, bottom.width, bottom.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, ptr(bottom.data));

	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

	glBindTexture(GL_TEXTURE_CUBE_MAP, 0);

	return cubemap;
}