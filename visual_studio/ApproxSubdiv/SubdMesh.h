#pragma once
#include "common.h"
#include "Geometry/TriangleMesh.h"

struct BufferTrait
{
	const void* data;
	int32_t size;
	int32_t offset;
	int32_t stride;
};

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

	void getPatch(BufferTrait &trait) const;
private:
	void initPatch();

private:
	static const uint32_t sPatchSize = 16;
	vector<Point3f> verts;
	vector<Point2f> uvs;
	vector<Normal3f> norms;
	vector<PolyIndex> fids;

	vector<Point3f> patch;
};

