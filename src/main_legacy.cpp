
#include "main.h"

int main()
{
	string filename;
#ifdef _DEBUG
	filename = "scene/obj/dragon2.obj";
#else
	cout << "Input filename: ";
	cin >> filename;
#endif // _DEBUG
	model_mesh = new SubdMesh(filename.c_str());

	GLFWwindow* window = nullptr;

	/* start GL context and O/S window using the GLFW helper library */
	if (!glfwInit()) {
		fprintf(stderr, "ERROR: could not start GLFW3\n");
		return 1;
	}

	window = glfwCreateWindow(win_width, win_height, "Carpenter", nullptr, nullptr);
	if (!window)
	{
		fprintf(stderr, "ERROR: could not open window with GLFW3\n");
		glfwTerminate();
		return 1;
	}
	glfwMakeContextCurrent(window);
	/* start GLEW extension handler */
	initGL();

	glfwSetKeyCallback(window, key_callback);
	//glfwSetWindowSize(window, win_width, win_height);
	glfwSetWindowSizeCallback(window, window_resize_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetScrollCallback(window, scroll_callback);

	GLFWcursor* cursor = glfwCreateStandardCursor(GLFW_ARROW_CURSOR);
	if (cursor) glfwSetCursor(window, cursor);
	glfwSetCursorEnterCallback(window, cursor_enter_callback);
	glfwSetCursorPosCallback(window, cursor_position_callback);

	while (!glfwWindowShouldClose(window))
	{
		if (view_changed)
		{
			window_refresh_callback(window);
			view_changed = false;
		}
		glfwPollEvents();
	}

	//glfwDestroyCursor(cursor);
	/* close GL context and any other GLFW resources */
	glfwTerminate();
	return 0;
}
