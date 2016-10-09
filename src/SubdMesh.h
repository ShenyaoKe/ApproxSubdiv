#pragma once
#include "common.h"
#include "Geometry/TriangleMesh.h"
#include "halfedge.h"

struct BufferTrait
{
	const void* data = nullptr;
	uint32_t count = 0;
	uint32_t size = 0;
	uint32_t offset = 0;
	uint32_t stride = 0;
};

struct BezierPatch
{
	Point3f patch[16];
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
	void savePatch() const;
private:
	void initPatch();

private:
	static const uint32_t sPatchSize = 16;
	vector<Point3f> verts;
	vector<Point2f> uvs;
	vector<Normal3f> norms;
	vector<PolyIndex> fids;

	unique_ptr<HDS_Mesh> mHDSMesh;
	vector<Point3f> patch4;
};

