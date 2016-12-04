#pragma once
#include "common.h"
#include "Geometry/TriangleMesh.h"
#include "halfedge.h"

struct BufferTrait
{
    const void* data    = nullptr;
    uint32_t    count   = 0;
    uint32_t    size    = 0;
    uint32_t    offset  = 0;
    uint32_t    stride  = 0;
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
	void getPatch(BufferTrait &bezier_trait,
					BufferTrait &gregory_trait) const;
	void savePatch() const;
private:
	Point3f faceCenter(uint32_t fid) const;
	Point3f edgeCenter(uint32_t heid) const;

	void initPatch();
	void genBezierPatch(
		const vector<bool> &bad_face_hash
	);
	void genGregoryPatch(
		vector<uint32_t> &vValenceCount,
		vector<uint32_t> &irreg_faces
	);
    void evalGregory() const;

private:
	static const uint32_t sBezierPatchSize = 16;
	static const uint32_t sGregoryPatchSize = 20;

	vector<Point3f> verts;
	vector<Point2f> uvs;
	vector<Normal3f> norms;
	vector<PolyIndex> fids;

	unique_ptr<HDS_Mesh> mHDSMesh;
	// Bezier Patch
	// 12---13---14---15
	// |    |    |    |
	// 08---09---10---11
	// |    |    |    |
	// 04---05---06---07
	// |    |    |    |
	// 00---01---02---03
	vector<Point3f> bezier_patch;
	// Gregory Patch Layout
	// 15-----17-----11-----10
	// |      |      |      |
	// |      19     13     |
	// 16--18          14--12
	// |                    |
	// |                    |
	// 02--04          08--06
	// |      03     09     |
	// |      |      |      |
	// 00-----01-----07-----05
	vector<Point3f> gregory_patch;
};

