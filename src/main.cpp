
#include "main.h"

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
	glfwSetScrollCallback(window, scroll_callback);

		

	GLFWcursor* cursor = glfwCreateStandardCursor(GLFW_ARROW_CURSOR);
	if (cursor) glfwSetCursor(window, cursor);
	glfwSetCursorPosCallback(window, cursor_position_callback);
	/* tell GL to only draw onto a pixel if the shape is closer to the viewer
	than anything already drawn at that pixel */
	

	// Setup buffers

	while (!glfwWindowShouldClose(window))
	{
		/* wipe the drawing surface clear */
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glClearColor(0.6, 0.6, 0.6, 1.0);

		if (!draw_wireframe)
		{
			//glEnable(GL_CULL_FACE);
			//glCullFace(GL_BACK); // cull back face
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
		glBindVertexArray(patch_vao);
		patch_shader->use_program();
		// Draw something
		glUniformMatrix4fv((*patch_shader)["view_matrix"], 1, GL_FALSE, view_cam->world_to_cam());
		glUniformMatrix4fv((*patch_shader)["proj_matrix"], 1, GL_FALSE, view_cam->cam_to_screen());
		glUniform1f((*patch_shader)["segments"], tess_seg);
		glPatchParameteri(GL_PATCH_VERTICES, 16);
		glDrawArrays(GL_PATCHES, 0, patchTrats.count);

		//glDisable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glBindVertexArray(model_vao);
		model_shader->use_program();
		// Draw something
		glUniformMatrix4fv((*model_shader)["view_matrix"], 1, GL_FALSE, view_cam->world_to_cam());
		glUniformMatrix4fv((*model_shader)["proj_matrix"], 1, GL_FALSE, view_cam->cam_to_screen());
		glDrawElements(GL_LINES_ADJACENCY, model_idx.size(), GL_UNSIGNED_INT, 0);


		glfwPollEvents();
		glfwSwapBuffers(window);
	}

	glfwDestroyCursor(cursor);
	/* close GL context and any other GLFW resources */
	glfwTerminate();
	return 0;
}
