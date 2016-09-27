#include <GL/glew.h> /* include GLEW and new version of GL on Windows */
#include <GL/glfw3.h> /* GLFW helper library */
#include "GLSLProgram.h"
#include "OGLViewer.h"
#include "SubdMesh.h"

static uint32_t win_width = 640;
static uint32_t win_height = 480;

unique_ptr<perspCamera> view_cam = make_unique<perspCamera>(
	Point3f(10, 6, 10), Point3f(0, 0, 0), Vector3f(0, 1, 0),
	win_width / static_cast<Float>(win_height));

SubdMesh model_mesh("scene/obj/torus_low.obj");
vector<GLfloat> model_verts;// vertices vbo
vector<GLuint> model_idx;// Normal coordinates vbo
BufferTrait patchTrats;
GLuint model_vert_vbo, model_ibo, model_vao;
unique_ptr<GLSLProgram> model_shader;

const GLubyte* renderer;
const GLubyte* version;

bool draw_wireframe = true;

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
	glDeleteBuffers(1, &model_vert_vbo);
	glDeleteBuffers(1, &model_ibo);
	glDeleteVertexArrays(1, &model_vao);

	glCreateBuffers(1, &model_vert_vbo);
	glNamedBufferData(model_vert_vbo, patchTrats.size, patchTrats.data, GL_STATIC_DRAW);

	// IBO
	/*glCreateBuffers(1, &model_ibo);
	glNamedBufferData(model_ibo, sizeof(GLuint) * model_idx.size(), &model_idx[0], GL_STATIC_DRAW);*/

	// VAO
	glCreateVertexArrays(1, &model_vao);
	glEnableVertexArrayAttrib(model_vao, 0);

	// Setup the formats
	glVertexArrayAttribFormat(model_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayVertexBuffer(model_vao, 0, model_vert_vbo, patchTrats.offset, patchTrats.stride);
	glVertexArrayAttribBinding(model_vao, 0, 0);
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

	model_mesh.exportIndexedVBO(&model_verts, nullptr, nullptr, &model_idx);
	bindMesh();

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_F && action == GLFW_PRESS)
	{
		view_cam.reset(new perspCamera(
			Point3f(10, 6, 10), Point3f(0, 0, 0), Vector3f(0, 1, 0),
			win_width / static_cast<Float>(win_height)));
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
	}
}

int main()
{
	GLFWwindow* window = nullptr;

	/* start GL context and O/S window using the GLFW helper library */
	if (!glfwInit()) {
		fprintf(stderr, "ERROR: could not start GLFW3\n");
		return 1;
	}

	window = glfwCreateWindow(win_width, win_height, "ApproxSubdiv", nullptr, nullptr);
	if (!window) {
		fprintf(stderr, "ERROR: could not open window with GLFW3\n");
		glfwTerminate();
		return 1;
	}
	glfwMakeContextCurrent(window);
	/* start GLEW extension handler */
	initGL();

	glfwSetKeyCallback(window, key_callback);
	glfwSetWindowSizeCallback(window, window_resize_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);

	/* tell GL to only draw onto a pixel if the shape is closer to the viewer
	than anything already drawn at that pixel */
	

	// Setup buffers

	while (!glfwWindowShouldClose(window))
	{
		/* wipe the drawing surface clear */
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glClearColor(0.6, 0.6, 0.6, 1.0);

		if (draw_wireframe)
		{
			glDisable(GL_CULL_FACE);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
		else
		{
			glEnable(GL_CULL_FACE);
			glCullFace(GL_BACK); // cull back face
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}

		glBindVertexArray(model_vao);
		model_shader->use_program();
		// Draw something
		glUniformMatrix4fv((*model_shader)["view_matrix"], 1, GL_FALSE, view_cam->world_to_cam());
		glUniformMatrix4fv((*model_shader)["proj_matrix"], 1, GL_FALSE, view_cam->cam_to_screen());
		glDrawElements(GL_LINES_ADJACENCY, model_idx.size(), GL_UNSIGNED_INT, 0);

		glfwPollEvents();
		glfwSwapBuffers(window);
	}

	/* close GL context and any other GLFW resources */
	glfwTerminate();
	return 0;
}
