#pragma once
#ifndef __OGLVIEWER__
#define __OGLVIEWER__

#include "GL/glew.h"
#include "common.h"

#include <QOpenGLWidget>
#include <QKeyEvent>
#include <QTimer>
#include <QTime>
#include <QString>
#include <QFileDialog>
#include <QOpenGLVertexArrayObject>
//#include <QGLFunctions>

#include "OpenGL_Utils/GLSLProgram.h"

#include "Math/MathUtil.h"
#include "Image/ImageData.h"
#include "Geometry/TriangleMesh.h"
//#include "Math/Matrix4D.h"
#include "Camera/perspCamera.h"
#include "SubdMesh.h"

//static int model_mat_loc;// Uniform matrix location
static GLfloat model_mat[16];// Uniform matrix buffer
//static int view_mat_loc;// Uniform matrix location
static GLfloat view_mat[16];
//static int proj_mat_loc;// Porjection matrix location
static GLfloat proj_mat[16];

const double eps = 5e-4;


// OpenGL Window in Qt
class OGLViewer : public QOpenGLWidget
{
	Q_OBJECT
public:
	enum Select_Mode
	{
		OBJECT_SELECT,
		COMPONENT_SELECT,
		FACE_COMPONENT_SELECT
	};
	OGLViewer(QWidget *parent = nullptr);
	~OGLViewer();

	//void update();
	public slots:
	void resetCamera();
	void initParas();
protected:
	void initializeGL() Q_DECL_OVERRIDE;
	void paintGL() Q_DECL_OVERRIDE;
	void resizeGL(int w, int h) Q_DECL_OVERRIDE;

	void keyPressEvent(QKeyEvent *e) Q_DECL_OVERRIDE;
	void mousePressEvent(QMouseEvent *e) Q_DECL_OVERRIDE;
	void mouseReleaseEvent(QMouseEvent *e) Q_DECL_OVERRIDE;
	void mouseMoveEvent(QMouseEvent *e) Q_DECL_OVERRIDE;
private:
	void bindMesh();
	void saveFrameBuffer();
public:
	double process_fps;
protected:
	unique_ptr<perspCamera> view_cam;
	Matrix4x4 proj, view;
private:
	int fps;
	int tcount;

	QTime process_time;
	int m_lastMousePos[2];
	int m_selectMode;
private: // OpenGL variables
	int display_mode = 0;

	SubdMesh model_mesh;
	vector<GLfloat> model_verts;// vertices vbo
	vector<GLuint> model_idx;// Normal coordinates vbo
	GLuint model_vert_vbo, model_ibo, model_vao;

	vector<GLfloat> filmgate, resgate;

	unique_ptr<GLSLProgram> model_shader;// OpenGL shader program

	friend class MainWindow;
};

#endif // __OGLVIEWER__