#pragma once

#include "GL/glew.h"
#include "common.h"

#include <QOpenGLWidget>
#include <QKeyEvent>
#include <QTimer>
#include <QTime>
#include <QString>
#include <QFileDialog>

#include "OpenGL_Utils/GLSLProgram.h"
#include "Camera/PerspectiveCamera.h"
#include "SubdMesh.h"

using namespace Kaguya;

class OglViewer : QOpenGLWidget
{
	Q_OBJECT
public:
	OglViewer(QWidget *parent = nullptr);
	~OglViewer();

	void resetCamera();
protected:
	void initializeGL() override;
	void paintGL() override;
	void resizeGL(int w, int h) override;

	void keyPressEvent(QKeyEvent *e) override;
	void mousePressEvent(QMouseEvent *e) override;
	void mouseReleaseEvent(QMouseEvent *e) override;
	void mouseMoveEvent(QMouseEvent *e) override;
private:
	void bindMesh();
	void bindPatchVertex();
	void bindBezierPatch();
	void bindQuadGregoryPatch();
	void bindTriGregoryPatch();

	GLuint createRenderObject(const RenderBufferTrait &trait);
	void saveFrameBuffer();
public:
protected:
private:
	float mProcessFps;
	int mFps;
	int mTimeCount;

	QTime mProcTime;
	int mLastMousePos[2];
	int mSelectMode;

	// OpenGL variables
	int mDisplayMode = 0;

	GLfloat mTessSeg = 4;
	bool mDrawCage = false;
	bool draw_wireframe = false;

	std::shared_ptr<SubdMesh> model_mesh;
	std::vector<GLfloat> model_verts;// vertices vbo
	std::vector<GLuint> model_idx;// Normal coordinates vbo
	GLuint model_vert_vbo, model_ibo, model_vao;
	std::unique_ptr<GLSLProgram> model_shader;

	GLuint mPatchVBO;
	BufferTrait mPatchVertexTrait;
	// Bezier Patch Buffer Objects
	GLuint mBezierPatchIBO, mBezierPatchVAO;
	BufferTrait mBezierPatchTrait;
	std::unique_ptr<GLSLProgram> mBezierPatchShader;
	// Gregory Patch Buffer Objects
	GLuint mQuadGregoryPatchIBO, mQuadGregoryPatchVAO;
	BufferTrait mQuadGregoryPatchTrait;
	std::unique_ptr<GLSLProgram> mQuadGregoryPatchShader;
	// Gregory Patch Buffer Objects
	GLuint mTriGregoryPatchIBO, mTriGregoryPatchVAO;
	BufferTrait mTriGregoryPatchTrait;
	std::unique_ptr<GLSLProgram> mTriGregoryPatchShader;

	vector<GLfloat> filmgate, resgate;

	unique_ptr<Kaguya::PerspectiveCamera> mViewCamera;

	friend class MainWindow;
};

