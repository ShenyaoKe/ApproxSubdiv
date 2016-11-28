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

SubdMesh* model_mesh;
vector<GLfloat> model_verts;// vertices vbo
vector<GLuint> model_idx;// Normal coordinates vbo
GLuint model_vert_vbo, model_ibo, model_vao;
unique_ptr<GLSLProgram> model_shader;
// Bezier Patch Buffer Objects
GLuint bPatch_vert_vbo, bPatch_vao;
BufferTrait bezier_patch_trait;
unique_ptr<GLSLProgram> bPatch_shader;
// Gregory Patch Buffer Objects
GLuint gPatch_vert_vbo, gPatch_vao;
BufferTrait gregory_patch_trait;
unique_ptr<GLSLProgram> gPatch_shader;

const GLubyte* renderer;
const GLubyte* version;

bool view_changed = true;
static int drawtime = 0;
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


void bindBezierPatch()
{
	glDeleteBuffers(1, &bPatch_vert_vbo);
	glDeleteVertexArrays(1, &bPatch_vao);

	glCreateBuffers(1, &bPatch_vert_vbo);
	glNamedBufferData(
		bPatch_vert_vbo,
		bezier_patch_trait.size,
		bezier_patch_trait.data,
		GL_STATIC_DRAW
	);

	// VAO
	glCreateVertexArrays(1, &bPatch_vao);
	glEnableVertexArrayAttrib(bPatch_vao, 0);

	// Setup the formats
	glVertexArrayAttribFormat(bPatch_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayVertexBuffer(
		bPatch_vao,
		0,
		bPatch_vert_vbo,
		bezier_patch_trait.offset,
		bezier_patch_trait.stride
	);
	glVertexArrayAttribBinding(bPatch_vao, 0, 0);
}

void bindGregoryPatch()
{
    glDeleteBuffers(1, &gPatch_vert_vbo);
    glDeleteVertexArrays(1, &gPatch_vao);

    glCreateBuffers(1, &gPatch_vert_vbo);
    glNamedBufferData(
        gPatch_vert_vbo,
        gregory_patch_trait.size,
        gregory_patch_trait.data,
        GL_STATIC_DRAW
    );

    // VAO
    glCreateVertexArrays(1, &gPatch_vao);
    glEnableVertexArrayAttrib(gPatch_vao, 0);

    // Setup the formats
    glVertexArrayAttribFormat(gPatch_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
    glVertexArrayVertexBuffer(
        gPatch_vao,
        0,
        gPatch_vert_vbo,
        gregory_patch_trait.offset,
        gregory_patch_trait.stride
    );
    glVertexArrayAttribBinding(gPatch_vao, 0, 0);
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
	bPatch_shader = make_unique<GLSLProgram>(
		"shaders/patch_vs.glsl",
		"shaders/patch_fs.glsl",
		nullptr,
		"shaders/patch_tc.glsl",
		"shaders/patch_te.glsl");
    gPatch_shader = make_unique<GLSLProgram>(
        "shaders/patch_vs.glsl",
        "shaders/patch_fs.glsl",
        "shaders/gregory_patch_gs.glsl",
        "shaders/gregory_patch_tc.glsl",
        "shaders/gregory_patch_te.glsl");

	model_mesh->exportIndexedVBO(&model_verts, nullptr, nullptr, &model_idx);
	bindMesh();
	model_mesh->getPatch(bezier_patch_trait, gregory_patch_trait);
	bindBezierPatch();
    bindGregoryPatch();

	glEnable(GL_BLEND);
	glEnable(GL_DEPTH_TEST); // enable depth-testing
	glDepthFunc(GL_LEQUAL);
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(0.6, 0.6, 0.6, 1.0);
}

void window_refresh_callback(GLFWwindow* window)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if (!draw_wireframe)
	{
		//glEnable(GL_CULL_FACE);
		//glCullFace(GL_BACK); // cull back face
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	else
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
    
	// Draw Bezier Patches
	glBindVertexArray(bPatch_vao);
	bPatch_shader->use_program();
	glUniformMatrix4fv((*bPatch_shader)["view_matrix"], 1, GL_FALSE, view_cam->world_to_cam());
	glUniformMatrix4fv((*bPatch_shader)["proj_matrix"], 1, GL_FALSE, view_cam->cam_to_screen());
	glUniform1f((*bPatch_shader)["segments"], tess_seg);
	glPatchParameteri(GL_PATCH_VERTICES, 16);
	glDrawArrays(GL_PATCHES, 0, bezier_patch_trait.count);

	if (draw_cage)
	{
		//glDisable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glBindVertexArray(model_vao);
		model_shader->use_program();
		// Draw something
		glUniformMatrix4fv((*model_shader)["view_matrix"], 1, GL_FALSE, view_cam->world_to_cam());
		glUniformMatrix4fv((*model_shader)["proj_matrix"], 1, GL_FALSE, view_cam->cam_to_screen());
		glDrawElements(GL_LINES_ADJACENCY, model_idx.size(), GL_UNSIGNED_INT, 0);
	}
    
    // Draw Gregory Patches
    glBindVertexArray(gPatch_vao);
    gPatch_shader->use_program();
    glUniformMatrix4fv((*gPatch_shader)["view_matrix"], 1, GL_FALSE, view_cam->world_to_cam());
    glUniformMatrix4fv((*gPatch_shader)["proj_matrix"], 1, GL_FALSE, view_cam->cam_to_screen());
    glUniform1f((*gPatch_shader)["segments"], tess_seg);
    glPatchParameteri(GL_PATCH_VERTICES, 20);
    glDrawArrays(GL_PATCHES, 0, gregory_patch_trait.count);

	glfwSwapBuffers(window);
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

	if (key == GLFW_KEY_Z && action == GLFW_PRESS)
	{
		zoom_camera = !zoom_camera;
	}

	if (key == GLFW_KEY_UP && action != GLFW_PRESS)
	{
		view_cam->zoom(0, 0.05, 0.0);
	}
	if (key == GLFW_KEY_DOWN && action != GLFW_RELEASE)
	{
		view_cam->zoom(0, -0.05, 0.0);
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
	view_changed = true;
}

void window_resize_callback(GLFWwindow* window, int width, int height)
{
	win_width = width;
	win_height = height;
	view_cam->resizeViewport(win_width / static_cast<Float>(win_height));
	glViewport(0, 0, win_width, win_height);

	view_changed = true;
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

	view_changed = true;
}

void cursor_enter_callback(GLFWwindow* window, int entered)
{
	if (entered)
	{
		glfwGetCursorPos(window, &cursor_last_x, &cursor_last_y);
	}
	else
	{
		// The cursor left the client area of the window
	}
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	double dx = xpos - cursor_last_x;
	double dy = ypos - cursor_last_y;

	if (move_camera)
	{
		view_cam->rotate(dy * 0.25f, -dx * 0.25f, 0.0);

		view_changed = true;
	}
	/*if (zoom_camera && dx != cursor_last_x)// zooming
	{
		view_cam->zoom(0.0, 0.0, dx * 0.05);

		view_changed = true;
	}
	if (!zoom_camera && dx != cursor_last_x && dy != cursor_last_y)
	{
		view_cam->zoom(-dx * 0.05, dy * 0.05, 0.0);
	}*/

	cursor_last_x = xpos;
	cursor_last_y = ypos;

}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	view_cam->zoom(0, 0, yoffset);

	view_changed = true;
}