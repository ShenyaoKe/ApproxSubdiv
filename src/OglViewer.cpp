#include "OglViewer.h"

OglViewer::OglViewer(QWidget* parent)
	: QOpenGLWidget(parent)
	, mViewCamera(new PerspectiveCamera(
		Kaguya::Point3f(6, 3, 6),
		Kaguya::Point3f(0, 1, 0),
		Kaguya::Vector3f(0, 1, 0),
		width() / float(height())))
	, model_mesh(new SubdMesh("scene/obj/pyramid.obj"))
{
}

OglViewer::~OglViewer()
{
}

/************************************************************************/
/* OpenGL Rendering Modules                                             */
/************************************************************************/
void OglViewer::initializeGL()
{
	// OpenGL extension initialization
	glewInit();
	// Print OpenGL version
	cout << "Renderer: " << glGetString(GL_RENDERER) << endl;
	cout << "OpenGL version supported " << glGetString(GL_VERSION) << endl;

	// Enable OpenGL features
	glEnable(GL_MULTISAMPLE);
	glEnable(GL_BLEND);
	glEnable(GL_DEPTH_TEST); // enable depth-testing
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glFrontFace(GL_CCW); // set counter-clock-wise vertex order to mean the front

	glClearColor(0.6, 0.6, 0.6, 0.0);

	// Create shader files
	model_shader = make_unique<GLSLProgram>("shaders/quad_vs.glsl",
											"shaders/quad_fs.glsl",
											"shaders/quad_gs.glsl");
	mBezierPatchShader = make_unique<GLSLProgram>("shaders/patch_vs.glsl",
												  "shaders/patch_fs.glsl",
												  nullptr,
												  "shaders/bezier_patch_tc.glsl",
												  "shaders/bezier_patch_te.glsl");
	mQuadGregoryPatchShader = make_unique<GLSLProgram>("shaders/patch_vs.glsl",
													   "shaders/patch_fs.glsl",
													   nullptr,
													   "shaders/gregory_patch_tc.glsl",
													   "shaders/gregory_patch_te.glsl");
	mTriGregoryPatchShader = make_unique<GLSLProgram>("shaders/patch_vs.glsl",
													  "shaders/patch_fs.glsl",
													  nullptr,
													  "shaders/tri_gregory_patch_tc.glsl",
													  "shaders/tri_gregory_patch_te.glsl");
	// Export vbo for shaders

	//model_mesh->exportIndexedVBO(&model_verts, nullptr, nullptr, &model_idx);
	//bindMesh();
	model_mesh->getPatch(mPatchVertexTrait,
						 mBezierPatchTrait,
						 mQuadGregoryPatchTrait,
						 mTriGregoryPatchTrait);
	bindPatchVertex();
	bindBezierPatch();
	bindQuadGregoryPatch();
	bindTriGregoryPatch();
}

void OglViewer::bindMesh()
{
	glDeleteBuffers(1, &model_vert_vbo);
	glDeleteBuffers(1, &model_ibo);
	glDeleteVertexArrays(1, &model_vao);

	glCreateBuffers(1, &model_vert_vbo);
	glNamedBufferData(model_vert_vbo,
					  sizeof(GLfloat) * model_verts.size(),
					  &model_verts[0],
					  GL_STATIC_DRAW);

	// IBO
	glCreateBuffers(1, &model_ibo);
	glNamedBufferData(model_ibo,
					  sizeof(GLuint) * model_idx.size(),
					  &model_idx[0],
					  GL_STATIC_DRAW);

	// VAO
	glCreateVertexArrays(1, &model_vao);
	glEnableVertexArrayAttrib(model_vao, 0);

	// Setup the formats
	glVertexArrayAttribFormat(model_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayVertexBuffer(model_vao, 0,
							  model_vert_vbo, 0,
							  sizeof(GLfloat) * 3);
	glVertexArrayAttribBinding(model_vao, 0, 0);

	glVertexArrayElementBuffer(model_vao, model_ibo);
}

void OglViewer::bindPatchVertex()
{
	glCreateBuffers(1, &mPatchVBO);
	glNamedBufferData(mPatchVBO,
					  mPatchVertexTrait.size,
					  mPatchVertexTrait.data,
					  GL_STATIC_DRAW);
}

void OglViewer::bindBezierPatch()
{
	glDeleteVertexArrays(1, &mBezierPatchVAO);

	glCreateBuffers(1, &mBezierPatchIBO);
	glNamedBufferData(mBezierPatchIBO,
					  mBezierPatchTrait.size,
					  mBezierPatchTrait.data,
					  GL_STATIC_DRAW);

	// VAO
	glCreateVertexArrays(1, &mBezierPatchVAO);
	glEnableVertexArrayAttrib(mBezierPatchVAO, 0);

	// Setup the formats
	glVertexArrayAttribFormat(mBezierPatchVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayVertexBuffer(mBezierPatchVAO, 0,
							  mPatchVBO,
							  mPatchVertexTrait.offset,
							  mPatchVertexTrait.stride);
	glVertexArrayAttribBinding(mBezierPatchVAO, 0, 0);
	glVertexArrayElementBuffer(mBezierPatchVAO, mBezierPatchIBO);
}

void OglViewer::bindQuadGregoryPatch()
{
	glDeleteBuffers(1, &mQuadGregoryPatchIBO);
	glDeleteVertexArrays(1, &mQuadGregoryPatchVAO);

	glCreateBuffers(1, &mQuadGregoryPatchIBO);
	glNamedBufferData(mQuadGregoryPatchIBO,
					  mQuadGregoryPatchTrait.size,
					  mQuadGregoryPatchTrait.data,
					  GL_STATIC_DRAW);

	// VAO
	glCreateVertexArrays(1, &mQuadGregoryPatchVAO);
	glEnableVertexArrayAttrib(mQuadGregoryPatchVAO, 0);

	// Setup the formats
	glVertexArrayAttribFormat(mQuadGregoryPatchVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayVertexBuffer(mQuadGregoryPatchVAO, 0,
							  mPatchVBO,
							  mPatchVertexTrait.offset,
							  mPatchVertexTrait.stride);
	glVertexArrayAttribBinding(mQuadGregoryPatchVAO, 0, 0);
	glVertexArrayElementBuffer(mQuadGregoryPatchVAO, mQuadGregoryPatchIBO);
}

void OglViewer::bindTriGregoryPatch()
{
	glDeleteBuffers(1, &mTriGregoryPatchIBO);
	glDeleteVertexArrays(1, &mTriGregoryPatchVAO);

	glCreateBuffers(1, &mTriGregoryPatchIBO);
	glNamedBufferData(mTriGregoryPatchIBO,
					  mTriGregoryPatchTrait.size,
					  mTriGregoryPatchTrait.data,
					  GL_STATIC_DRAW);

	// VAO
	glCreateVertexArrays(1, &mTriGregoryPatchVAO);
	glEnableVertexArrayAttrib(mTriGregoryPatchVAO, 0);

	// Setup the formats
	glVertexArrayAttribFormat(mTriGregoryPatchVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayVertexBuffer(mTriGregoryPatchVAO, 0,
							  mPatchVBO,
							  mPatchVertexTrait.offset,
							  mPatchVertexTrait.stride);
	glVertexArrayAttribBinding(mTriGregoryPatchVAO, 0, 0);
	glVertexArrayElementBuffer(mTriGregoryPatchVAO, mTriGregoryPatchIBO);
}

GLuint OglViewer::createRenderObject(const RenderBufferTrait &trait)
{
	GLuint vbo, ibo, vao;

	// VBO
	glCreateBuffers(1, &vbo);
	glNamedBufferData(vbo, trait.vertex.size, trait.vertex.data, GL_STATIC_DRAW);
	// IBO
	glCreateBuffers(1, &ibo);
	glNamedBufferData(ibo, trait.index.size, trait.index.data, GL_STATIC_DRAW);
	//indexCount = trait.index.count;

	// Bind VAO
	glCreateVertexArrays(1, &vao);
	glEnableVertexArrayAttrib(vao, 0);

	// Attach VBO and IBO to VAO
	glVertexArrayAttribFormat(vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayVertexBuffer(vao, 0, vbo, trait.vertex.offset, trait.vertex.stride);
	glVertexArrayAttribBinding(vao, 0, 0);
	glVertexArrayElementBuffer(vao, ibo);

	return vao;
}

void OglViewer::paintGL()
{
	// Make current window
	makeCurrent();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_BACK); // cull back face
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

	if (mDrawCage)
	{
		//glDisable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glBindVertexArray(model_vao);
		model_shader->use_program();
		// Draw something
		glUniformMatrix4fv((*model_shader)["view_matrix"], 1, GL_FALSE,
						   mViewCamera->world_to_cam());
		glUniformMatrix4fv((*model_shader)["proj_matrix"], 1, GL_FALSE,
						   mViewCamera->cam_to_screen());
		glDrawElements(GL_LINES_ADJACENCY, model_idx.size(), GL_UNSIGNED_INT, 0);
	}

	// Draw Bezier Patches
	glBindVertexArray(mBezierPatchVAO);
	mBezierPatchShader->use_program();
	glUniformMatrix4fv((*mBezierPatchShader)["view_matrix"], 1, GL_FALSE,
					   mViewCamera->world_to_cam());
	glUniformMatrix4fv((*mBezierPatchShader)["proj_matrix"], 1, GL_FALSE,
					   mViewCamera->cam_to_screen());
	glUniform1f((*mBezierPatchShader)["segments"], mTessSeg);
	glPatchParameteri(GL_PATCH_VERTICES, SubdMesh::sQuadBezierPatchSize);
	glDrawElements(GL_PATCHES, mBezierPatchTrait.count, GL_UNSIGNED_INT, 0);


	// Draw Quadrilateral Gregory Patches
	glBindVertexArray(mQuadGregoryPatchVAO);
	mQuadGregoryPatchShader->use_program();
	glUniformMatrix4fv((*mQuadGregoryPatchShader)["view_matrix"], 1, GL_FALSE,
					   mViewCamera->world_to_cam());
	glUniformMatrix4fv((*mQuadGregoryPatchShader)["proj_matrix"], 1, GL_FALSE,
					   mViewCamera->cam_to_screen());
	glUniform1f((*mQuadGregoryPatchShader)["segments"], mTessSeg);
	glPatchParameteri(GL_PATCH_VERTICES, SubdMesh::sQuadGregoryPatchSize);
	glDrawElements(GL_PATCHES, mQuadGregoryPatchTrait.count, GL_UNSIGNED_INT, 0);

	// Draw Triangular Gregory Patches
	glBindVertexArray(mTriGregoryPatchVAO);
	mTriGregoryPatchShader->use_program();
	glUniformMatrix4fv((*mTriGregoryPatchShader)["view_matrix"], 1, GL_FALSE,
					   mViewCamera->world_to_cam());
	glUniformMatrix4fv((*mTriGregoryPatchShader)["proj_matrix"], 1, GL_FALSE,
					   mViewCamera->cam_to_screen());
	glUniform1f((*mTriGregoryPatchShader)["segments"], mTessSeg);
	glPatchParameteri(GL_PATCH_VERTICES, SubdMesh::sTriGregoryPatchSize);
	glDrawElements(GL_PATCHES, mTriGregoryPatchTrait.count, GL_UNSIGNED_INT, 0);

	glBindVertexArray(0);
}
// Redraw function

// Resize function
void OglViewer::resizeGL(int /*w*/, int /*h*/)
{
	// Widget resize operations
	mViewCamera->resizeViewport(width() / Float(height()));
}
/************************************************************************/
/* Qt User Operation Functions                                          */
/************************************************************************/
void OglViewer::keyPressEvent(QKeyEvent *e)
{
	if (e->key() == Qt::Key_W)
	{
		mViewCamera->zoom(0.0, 0.0, 0.10);
	}
	else if (e->key() == Qt::Key_S)
	{
		mViewCamera->zoom(0.0, 0.0, -0.10);
	}
	else if (e->key() == Qt::Key_Q)
	{
		mViewCamera->zoom(0.10, 0.0, 0.0);
	}
	else if (e->key() == Qt::Key_A)
	{
		mViewCamera->zoom(-0.10, 0.0, 0.0);
	}
	else if (e->key() == Qt::Key_E)
	{
		mViewCamera->zoom(0.0, 0.10, 0.0);
	}
	else if (e->key() == Qt::Key_D)
	{
		mViewCamera->zoom(0.0, -0.10, 0.0);
	}
	else if (e->key() == Qt::Key_Plus)
	{
		mTessSeg++;
	}
	else if (e->key() == Qt::Key_Minus)
	{
		mTessSeg = std::max(mTessSeg - 1, 1.0f);
	}
	// Save frame buffer
	else if (e->key() == Qt::Key_P && e->modifiers() == Qt::ControlModifier)
	{
		saveFrameBuffer();
	}
	//////////////////////////////////////////////////////////////////////////
	else
	{
		QOpenGLWidget::keyPressEvent(e);
	}
	update();
}

void OglViewer::mousePressEvent(QMouseEvent *e)
{
	mLastMousePos[0] = e->x();
	mLastMousePos[1] = e->y();
	if ((e->buttons() == Qt::LeftButton) && (e->modifiers() == Qt::AltModifier))
	{
		// Do something here
	}
};

void OglViewer::mouseReleaseEvent(QMouseEvent *e)
{
	mLastMousePos[0] = e->x();
	mLastMousePos[1] = e->y();
}

void OglViewer::mouseMoveEvent(QMouseEvent *e)
{
	int dx = e->x() - mLastMousePos[0];
	int dy = e->y() - mLastMousePos[1];

	//printf("dx: %d, dy: %d\n", dx, dy);

	if ((e->buttons() == Qt::LeftButton) && (e->modifiers() == Qt::AltModifier))
	{
		mViewCamera->rotate(dy * 0.25, -dx * 0.25, 0.0);
		update();
	}
	else if ((e->buttons() == Qt::RightButton) && (e->modifiers() == Qt::AltModifier))
	{
		if (dx != e->x() && dy != e->y())
		{
			mViewCamera->zoom(0.0, 0.0, dx * 0.05);
			update();
		}
	}
	else if ((e->buttons() == Qt::MidButton) && (e->modifiers() == Qt::AltModifier))
	{
		if (dx != e->x() && dy != e->y())
		{
			mViewCamera->zoom(-dx * 0.05, dy * 0.05, 0.0);
			update();
		}
	}
	else
	{
		QOpenGLWidget::mouseMoveEvent(e);
	}

	mLastMousePos[0] = e->x();
	mLastMousePos[1] = e->y();
}
/************************************************************************/
/* Application Functions                                                */
/************************************************************************/
void OglViewer::resetCamera()
{
	mViewCamera.reset(new PerspectiveCamera(Point3f(10, 6, 10),
											Point3f(0, 0, 0),
											Vector3f(0, 1, 0),
											width() / Float(height())));
}

void OglViewer::saveFrameBuffer()
{
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Screenshot file...", "default", tr("PNG(*.png)"));
	grab().save(filename);
}
