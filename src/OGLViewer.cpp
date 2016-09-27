#include "OGLViewer.h"
#if 0

OGLViewer::OGLViewer(QWidget *parent)
	: m_selectMode(OBJECT_SELECT)
	, view_cam(new perspCamera(
		Point3f(10, 6, 10), Point3f(0, 0, 0), Vector3f(0, 1, 0),
		width() / static_cast<Float>(height())))
	, model_mesh("../../scene/obj/torus_low.obj")
{
	auto ctx = this->context();
	cout << "default fbo:\t" << defaultFramebufferObject() << endl;
	// Set surface format for current widget
	QSurfaceFormat format;
	format.setDepthBufferSize(32);
	format.setStencilBufferSize(8);
	//format.setSamples(4);
	format.setVersion(4, 5);
	format.setProfile(QSurfaceFormat::CoreProfile);
	this->setFormat(format);

	view_cam->exportVBO(view_mat, proj_mat, nullptr);
}

OGLViewer::~OGLViewer()
{
}
/************************************************************************/
/* OpenGL Rendering Modules                                             */
/************************************************************************/
void OGLViewer::initializeGL()
{
	// OpenGL extention initialization
	glewInit();
	// Print OpenGL vertion
	cout << "Renderer: " << glGetString(GL_RENDERER) << endl;
	cout << "OpenGL version supported " << glGetString(GL_VERSION) << endl;

	// Enable OpenGL features
	glEnable(GL_MULTISAMPLE);
	//glEnable(GL_LINE_SMOOTH);
	//glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_BLEND);
	glEnable(GL_DEPTH_TEST); // enable depth-testing
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glFrontFace(GL_CCW); // set counter-clock-wise vertex order to mean the front


	//////////////////////////////////////////////////////////////////////////

	// Create shader files
	model_shader.reset(new GLSLProgram(
		"shaders/quad_vs.glsl", "shaders/quad_fs.glsl",
		"shaders/quad_gs.glsl"));

	// Export vbo for shaders
	model_mesh.exportIndexedVBO(&model_verts, nullptr, nullptr, &model_idx);

	bindMesh();
}

void OGLViewer::bindMesh()
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

void OGLViewer::bindPatch()
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

	//glVertexArrayElementBuffer(model_vao, model_ibo);
}

void OGLViewer::paintGL()
{
	// Make curent window
	makeCurrent();

	glDisable(GL_MULTISAMPLE);
	// Clear background and color buffer
	glClearColor(0.6, 0.6, 0.6, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//////////////////////////////////////////////////////////////////////////
	// Model
	/*glDisable(GL_CULL_FACE);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);*/
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK); // cull back face
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glBindVertexArray(model_vao);
	model_shader->use_program();

	// Apply uniform matrix
	//glUniformMatrix4fv(model_mat_loc, 1, GL_FALSE, model_mat);
	glUniformMatrix4fv((*model_shader)["view_matrix"], 1, GL_FALSE, view_mat);
	glUniformMatrix4fv((*model_shader)["proj_matrix"], 1, GL_FALSE, proj_mat);
	glDrawElements(GL_LINES_ADJACENCY, model_idx.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);


//////////////////////////////////////////////////////////////////////////
	glPatchParameteri(GL_PATCH_VERTICES, 16);
	glDrawArrays(GL_PATCHES, hmsh_vtx_offset[i], patchSize);
}

// Resize function
void OGLViewer::resizeGL(int w, int h)
{
	// Widget resize operations
	view_cam->resizeViewport(width() / static_cast<double>(height()));
	view_cam->exportVBO(nullptr, proj_mat, nullptr);
}
/************************************************************************/
/* Qt User Operation Functions                                          */
/************************************************************************/
void OGLViewer::keyPressEvent(QKeyEvent *e)
{
	if (e->key() == Qt::Key_W)
	{
		view_cam->zoom(0.0, 0.0, 0.10);
		view_cam->exportVBO(view_mat, nullptr, nullptr);

	}
	else if (e->key() == Qt::Key_S)
	{
		view_cam->zoom(0.0, 0.0, -0.10);
		view_cam->exportVBO(view_mat, nullptr, nullptr);
	}
	else if (e->key() == Qt::Key_Q)
	{
		view_cam->zoom(0.10, 0.0, 0.0);
		view_cam->exportVBO(view_mat, nullptr, nullptr);
	}
	else if (e->key() == Qt::Key_A)
	{
		view_cam->zoom(-0.10, 0.0, 0.0);
		view_cam->exportVBO(view_mat, nullptr, nullptr);
	}
	else if (e->key() == Qt::Key_E)
	{
		view_cam->zoom(0.0, 0.10, 0.0);
		view_cam->exportVBO(view_mat, nullptr, nullptr);
	}
	else if (e->key() == Qt::Key_D)
	{
		view_cam->zoom(0.0, -0.10, 0.0);
		view_cam->exportVBO(view_mat, nullptr, nullptr);
	}
	else if (e->key() == Qt::Key_Home)
	{
		initParas();
	}
	// Save frame buffer
	else if (e->key() == Qt::Key_P && e->modifiers() == Qt::ControlModifier)
	{
		this->saveFrameBuffer();
	}
	//////////////////////////////////////////////////////////////////////////
	else
	{
		QOpenGLWidget::keyPressEvent(e);
	}
	update();
}

void OGLViewer::mousePressEvent(QMouseEvent *e)
{
	m_lastMousePos[0] = e->x();
	m_lastMousePos[1] = e->y();
	if ((e->buttons() == Qt::LeftButton) && (e->modifiers() == Qt::AltModifier))
	{
		// Do something here
	}
};

void OGLViewer::mouseReleaseEvent(QMouseEvent *e)
{
	m_lastMousePos[0] = e->x();
	m_lastMousePos[1] = e->y();
}

void OGLViewer::mouseMoveEvent(QMouseEvent *e)
{
	int dx = e->x() - m_lastMousePos[0];
	int dy = e->y() - m_lastMousePos[1];

	//printf("dx: %d, dy: %d\n", dx, dy);

	if ((e->buttons() == Qt::LeftButton) && (e->modifiers() == Qt::AltModifier))
	{
		view_cam->rotate(dy * 0.25, -dx * 0.25, 0.0);
		view_cam->exportVBO(view_mat, nullptr, nullptr);
		update();
	}
	else if ((e->buttons() == Qt::RightButton) && (e->modifiers() == Qt::AltModifier))
	{
		if (dx != e->x() && dy != e->y())
		{
			view_cam->zoom(0.0, 0.0, dx * 0.05);
			view_cam->exportVBO(view_mat, nullptr, nullptr);
			update();
		}
	}
	else if ((e->buttons() == Qt::MidButton) && (e->modifiers() == Qt::AltModifier))
	{
		if (dx != e->x() && dy != e->y())
		{
			view_cam->zoom(-dx * 0.05, dy * 0.05, 0.0);
			view_cam->exportVBO(view_mat, nullptr, nullptr);
			update();
		}
	}
	else
	{
		QOpenGLWidget::mouseMoveEvent(e);
	}

	m_lastMousePos[0] = e->x();
	m_lastMousePos[1] = e->y();
}
/************************************************************************/
/* Application Functions                                                */
/************************************************************************/
void OGLViewer::resetCamera()
{
	view_cam.reset(new perspCamera(
		Point3f(10, 6, 10), Point3f(0, 0, 0), Vector3f(0, 1, 0),
		width() / static_cast<Float>(height())));
	view_cam->exportVBO(view_mat, proj_mat, nullptr);
	//update();
}
void OGLViewer::initParas()
{
	update();
}

void OGLViewer::saveFrameBuffer()
{
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Screenshot file...", "default", tr("PNG(*.png)"));
	this->grab().save(filename);
}
#endif