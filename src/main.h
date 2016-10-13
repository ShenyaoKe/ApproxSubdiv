#pragma once
#include <GL/glew.h> /* include GLEW and new version of GL on Windows */
#include <GL/glfw3.h> /* GLFW helper library */
#include "GLSLProgram.h"
#include "OGLViewer.h"
#include "SubdMesh.h"


static int32_t win_width = 640;
static int32_t win_height = 480;
static double cursor_last_x, cursor_last_y;
static GLfloat tess_seg = 4;

unique_ptr<perspCamera> view_cam = make_unique<perspCamera>(
	Point3f(5, 3, 5), Point3f(0, 0, 0), Vector3f(0, 1, 0),
	win_width / static_cast<Float>(win_height));

SubdMesh model_mesh("scene/obj/dragon_scaled.obj");
vector<GLfloat> model_verts;// vertices vbo
vector<GLuint> model_idx;// Normal coordinates vbo
BufferTrait patchTrats;
GLuint model_vert_vbo, model_ibo, model_vao;
GLuint patch_vert_vbo, patch_vao;
unique_ptr<GLSLProgram> model_shader;
unique_ptr<GLSLProgram> patch_shader;

const GLubyte* renderer;
const GLubyte* version;

bool draw_cage = false;
bool draw_wireframe = false;
bool move_camera = false;
bool zoom_camera = false;

void bindMesh()
{
	glDeleteBuffers(1, &model_vert_vbo);
	glDeleteBuffers(1, &model_ibo);
	glDeleteVertexArrays(1, &model_vao);

	glCreateBuffers(1, &model_vert_vbo);
	glNamedBufferData(model_vert_vbo, sizeof(GLfloat) * model_verts.size(), &model_verts[0], GL_STATIC_DRAW);

	// IBO
	glCreateBuffers(1, &model_ibo);
	glNamedBufferData(model_ibo, sizeof(GLuint) * model_idx.size(), &model_idx[0], GL_STATIC_DRAW);

	// VAO
	glCreateVertexArrays(1, &model_vao);
	glEnableVertexArrayAttrib(model_vao, 0);

	// Setup the formats
	glVertexArrayAttribFormat(model_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayVertexBuffer(model_vao, 0, model_vert_vbo, 0, sizeof(GLfloat) * 3);
	glVertexArrayAttribBinding(model_vao, 0, 0);

	glVertexArrayElementBuffer(model_vao, model_ibo);
}


void bindPatch()
{
	glDeleteBuffers(1, &patch_vert_vbo);
	glDeleteVertexArrays(1, &patch_vao);

	glCreateBuffers(1, &patch_vert_vbo);
	glNamedBufferData(
		patch_vert_vbo,
		patchTrats.size,
		patchTrats.data,
		GL_STATIC_DRAW
	);

	// VAO
	glCreateVertexArrays(1, &patch_vao);
	glEnableVertexArrayAttrib(patch_vao, 0);

	// Setup the formats
	glVertexArrayAttribFormat(patch_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayVertexBuffer(
		patch_vao,
		0,
		patch_vert_vbo,
		patchTrats.offset,
		patchTrats.stride
	);
	glVertexArrayAttribBinding(patch_vao, 0, 0);
}

void initGL()
{
	renderer = glGetString(GL_RENDERER); /* get renderer string */
	version = glGetString(GL_VERSION); /* version as a string */
	printf("Renderer: %s\n", renderer);
	printf("OpenGL version supported %s\n", version);

	glewExperimental = GL_TRUE;
	glewInit();

	/* geometry to use. these are 3 xyz points (9 floats total) to make a triangle */
	model_shader = make_unique<GLSLProgram>(
		"shaders/quad_vs.glsl", "shaders/quad_fs.glsl", "shaders/quad_gs.glsl");
	patch_shader = make_unique<GLSLProgram>(
		"shaders/patch_vs.glsl",
		"shaders/patch_fs.glsl",
		nullptr,
		"shaders/patch_tc.glsl",
		"shaders/patch_te.glsl");

	model_mesh.exportIndexedVBO(&model_verts, nullptr, nullptr, &model_idx);
	bindMesh();
	model_mesh.getPatch(patchTrats);
	bindPatch();

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_F && action == GLFW_PRESS)
	{
		view_cam.reset(new perspCamera(
			Point3f(10, 6, 10), Point3f(0, 0, 0), Vector3f(0, 1, 0),
			win_width / static_cast<Float>(win_height)));
	}
	if (key == GLFW_KEY_W && action == GLFW_PRESS)
	{
		draw_wireframe = !draw_wireframe;
	}
	if (key == GLFW_KEY_C && action == GLFW_PRESS)
	{
		draw_cage = !draw_cage;
	}
	if (key == GLFW_KEY_KP_ADD)
	{
		tess_seg += 1.0f;
	}
	if (key == GLFW_KEY_KP_SUBTRACT)
	{
		tess_seg -= 1.0f;
		if (tess_seg < 1.0f)
		{
			tess_seg = 1.0f;
		}
	}
}

void window_resize_callback(GLFWwindow* window, int width, int height)
{
	win_width = width;
	win_height = height;
	view_cam->resizeViewport(win_width / static_cast<Float>(win_height));
	glViewport(0, 0, win_width, win_height);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
	{
		// Do sth
		zoom_camera = true;
	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
	{
		glfwGetCursorPos(window, &cursor_last_x, &cursor_last_y);
		move_camera = true;
	}
	if (action == GLFW_RELEASE)
	{
		move_camera = zoom_camera = false;
	}
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	double dx = xpos - cursor_last_x;
	double dy = ypos - cursor_last_y;

	if (move_camera)
	{
		view_cam->rotate(dy * 0.25f, -dx * 0.25f, 0.0);
	}
	if (zoom_camera && dx != cursor_last_x)// zooming
	{
		view_cam->zoom(0.0, 0.0, dx * 0.05);
	}
	/*if (!zoom_camera && dx != cursor_last_x && dy != cursor_last_y)
	{
		view_cam->zoom(-dx * 0.05, dy * 0.05, 0.0);
	}*/

	cursor_last_x = xpos;
	cursor_last_y = ypos;
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	view_cam->zoom(0, 0, yoffset);
}