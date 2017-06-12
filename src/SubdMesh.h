#pragma once
#include "common.h"
#include "Geometry/TriangleMesh.h"
#include "HalfEdge.h"

namespace Kaguya
{

struct BufferTrait
{
	const void* data = nullptr;
	uint32_t    count = 0;
	uint32_t    size = 0;
	uint32_t    offset = 0;
	uint32_t    stride = 0;
};

struct BezierPatch
{
	Point3f patch[16];
};
struct GregoryPatch
{
	Point3f patch[20];
};

struct GeomContext
{
	vector<uint32_t> mVertexValence;
	unordered_set<uint32_t> mBoundaryVertex, mCornerVertex;
	vector<bool> mBadFaceHash;

	vector<Point3f> mEdgeMiddle;
	vector<Point3f> mFaceCenter;
};

class SubdMesh
{
public:
	static const uint32_t sTriGregoryPatchSize = 15;
	static const uint32_t sQuadBezierPatchSize = 16;
	static const uint32_t sQuadGregoryPatchSize = 20;

	SubdMesh(const std::string &filename);
	~SubdMesh();

	void exportIndexedVBO(vector<Float>* vtx_array = nullptr,
						  vector<Float>* uv_array = nullptr,
						  vector<Float>* norm_array = nullptr,
						  vector<uint32_t>* idx_array = nullptr) const;

	void getPatch(BufferTrait &vertexBuffer,
				  BufferTrait &bezierPatchIndexBuffer,
				  BufferTrait &quadGregoryPatchIndexBuffer,
				  BufferTrait &triGregoryPatchIndexBuffer) const;
	void savePatch() const;

private:
	Point3f faceCenter(SizeType fid) const;
	Point3f edgeCenter(SizeType heid) const;

	void process();
	void specifyBoundaries(GeomContext &context);
	void specifyPatchType(GeomContext &context);
	// Cache edge middle points and face centers
	void cacheGeometry(GeomContext &context);

	// Generate Patch
	void computeCornerPoints(const GeomContext &context);
	void computeEdgePoints(const GeomContext &context);

	// Generate Patch topology and compute face points
	void computePatches(const GeomContext &context);
	void computeBezierPatch(SizeType fid);
	void computeQuadGregoryPatch(SizeType fid, const GeomContext &context);
	void computeTriGregoryPatch(SizeType fid, const GeomContext &context);

	Point3f computeQuadFacePoint(const Point3f &p0,
								 const Point3f &m_ip1, const Point3f &m_im1,
								 const Point3f &fc_i, const Point3f &fc_im1,
								 const Point3f &edge_near, const Point3f &edge_far,
								 SizeType valenceNear, SizeType valenceFar) const;
	Point3f computeBoundaryQuadFacePoint(const Point3f &p0,
										 const Point3f &nearP1,
										 const Point3f &nearP2,
										 const Point3f &diagP,
										 SizeType valence);

	void evalGregory() const;

#ifdef _DEBUG
	void hack();
#endif // _DEBUG


private:

	vector<Point3f> mVerts;
	vector<Point2f> uvs;
	vector<Normal3f> norms;
	vector<SizeType> mFaceIndices;
	vector<SizeType> mFaceSideCount;
	vector<SizeType> mFaceIdOffset;

	unique_ptr<HDS::Mesh> mHDSMesh;

	vector<Point3f> mPatchVertexBuffer;
	// Bezier Patch
	// 12---13---14---15
	// |    |    |    |
	// 08---09---10---11
	// |    |    |    |
	// 04---05---06---07
	// |    |    |    |
	// 00---01---02---03
	//vector<BezierPatch> mBezierPatches;
	vector<SizeType> mBezierPatchIndices;

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
	//vector<GregoryPatch> mGregoryPatch;
	vector<SizeType> mQuadGregoryPatchIndices;

	// Triangular Gregory Patch Layout
	//               10
	//             /    \
	//           11      12
	//          /  \    /  \
	//        /    13  14    \
	//      02--04        08--06
	//     /      03    09      \
	//   /         |    |         \
	//  00--------01----07--------05
	vector<SizeType> mTriGregoryPatchIndices;
};

}
