#pragma once

#include "vcl/vcl.hpp"
class WaterFrameBuffers {

private:
	static GLuint const REFLECTION_HEIGHT = 180;
	static GLuint const REFRACTION_HEIGHT = 720;

	GLuint reflectionFrameBuffer;
	GLuint reflectionTexture;
	GLuint reflectionDepthBuffer;

	GLuint refractionFrameBuffer;
	GLuint refractionTexture;
	GLuint refractionDepthTexture;

protected:
	static  GLuint const REFLECTION_WIDTH = 320;
	static  GLuint const REFRACTION_WIDTH = 1280;


private:
	void initialiseReflectionFrameBuffer(GLFWwindow* window);

	void initialiseRefractionFrameBuffer(GLFWwindow* window);

	void bindFrameBuffer(GLuint frameBuffer, GLuint width, GLuint height);
	GLuint createFrameBuffer();
	GLuint createTextureAttachment(GLuint width, GLuint height);
	GLuint createDepthTextureAttachment(GLuint width, GLuint height);
	GLuint createDepthBufferAttachment(GLuint width, GLuint height);


public:
	WaterFrameBuffers(GLFWwindow* window);

	void cleanUp();
	void bindReflectionFrameBuffer();

	void bindRefractionFrameBuffer();

	void unbindCurrentFrameBuffer(GLFWwindow* window);

	GLuint getReflectionTexture();

	GLuint getRefractionTexture();

	GLuint getRefractionDepthTexture();




};