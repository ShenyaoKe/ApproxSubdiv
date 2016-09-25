#pragma once
#include "common.h"
#include "Geometry/TriangleMesh.h"

class SubdMesh
{
public:
	SubdMesh(const char* filename);
	~SubdMesh();

	void exportIndexedVBO(
		vector<Float>* vtx_array = nullptr,
		vector<Float>* uv_array = nullptr,
		vector<Float>* norm_array = nullptr,
		vector<unsigned int>* idx_array = nullptr) const;

private:
	vector<Point3f> verts;
	vector<Point2f> uvs;
	vector<Normal3f> norms;
	vector<PolyIndex> fids;
};

