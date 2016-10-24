#pragma once
#include "common.h"
#include "Geometry/TriangleMesh.h"

static const int32_t sInvalidHDS = -1;

/*
* Vertex
*/
class HDS_Vertex
{
public:
	HDS_Vertex() : index(uid++), heid(sInvalidHDS) {}
	~HDS_Vertex() {}

	static void resetIndex() { uid = 0; }

	int32_t index;
	int32_t heid;
private:
	static int32_t uid;
};

/*
* Half-Edge
*/
class HDS_HalfEdge
{
public:
	static void resetIndex() { uid = 0; }
	static void matchIndexToSize(size_t	size) { uid = size; }

	HDS_HalfEdge() : index(uid++), fid(sInvalidHDS), vid(sInvalidHDS)
		, prev_offset(0), next_offset(0), flip_offset(0), isBoundary(false) {}
	~HDS_HalfEdge(){}

	//HDS_HalfEdge(const HDS_HalfEdge& other);
	//HDS_HalfEdge operator=(const HDS_HalfEdge& other);

	// Get the explicit pointer to corresponding edges
	HDS_HalfEdge* prev() { return this + prev_offset; }
	HDS_HalfEdge* next() { return this + next_offset; }
	HDS_HalfEdge* flip() { return this + flip_offset; }
	HDS_HalfEdge* rotCW() { return flip()->next(); }
	HDS_HalfEdge* rotCCW() { return prev()->flip(); }
	const HDS_HalfEdge* prev() const { return this + prev_offset; }
	const HDS_HalfEdge* next() const { return this + next_offset; }
	const HDS_HalfEdge* flip() const { return this + flip_offset; }
	const HDS_HalfEdge* rotCW() const { return flip()->next(); }
	const HDS_HalfEdge* rotCCW() const { return prev()->flip(); }

	void setFlip(HDS_HalfEdge* f_e)
	{ flip_offset = f_e - this; f_e->flip_offset = -flip_offset; }

	//////////////////////////////////////////////////////////////////////////
	int32_t index;
	int32_t fid;
	int32_t vid;

	// Offset to index of previous/nex/flip edge
	// previous/nex/flip edge doesn't exist
	// when (previous/nex/flip == 0)
	int32_t prev_offset, next_offset, flip_offset;

	bool isBoundary;
private:
	static int32_t uid;
};

/*
* Face
*/
class HDS_Face
{
public:
	static void resetIndex() { uid = 0; }

	HDS_Face() : index(uid++), heid(sInvalidHDS), isNullFace(false) {}
	~HDS_Face() {}

	// Get the connected half-edge id
	// Explicit pointer access is handled by HDS_Mesh
	int32_t heID() const { return heid; }

	// Member data
	int32_t index;
	int32_t heid;
	bool isNullFace;
private:
	static int32_t uid;
};

using vert_t = HDS_Vertex;
using he_t = HDS_HalfEdge;
using face_t = HDS_Face;
/*
* Mesh
*/
class HDS_Mesh
{
public:
	HDS_Mesh() = delete;
	HDS_Mesh(vector<vert_t> &vs, vector<he_t> &hes, vector<face_t> &fs)
		: verts(std::move(vs))
		, halfedges(std::move(hes))
		, faces(std::move(fs)) {}
	HDS_Mesh(const HDS_Mesh &other)
		: verts(other.verts), halfedges(other.halfedges), faces(other.faces) {}
	~HDS_Mesh() {}

	// Reset UID in each component
	// Mask(Bitwise Operation): face|edge|vertex
	//    e.g. All(111==3), Vertex Only(001==1), Vertex+Edge(011==3)
	static void resetIndex(uint8_t reset_mask = 7) {
		if (reset_mask & 1) HDS_Vertex::resetIndex();
		if (reset_mask & 2) HDS_HalfEdge::resetIndex();
		if (reset_mask & 4) HDS_Face::resetIndex();
	}
	
	void printInfo(const string &msg = "") {
		if (!msg.empty()) cout << msg << endl;
		cout << "#vertices = " << verts.size() << endl;
		cout << "#faces = " << faces.size() << endl;
		cout << "#half edges = " << halfedges.size() << endl;
	}

	he_t* heFromFace(int32_t fid) { return &halfedges[faces[fid].heid]; }
	he_t* heFromVert(int32_t vid) { return &halfedges[verts[vid].heid]; }
	vert_t* vertFromHe(int32_t heid) { return &verts[halfedges[heid].vid]; }
	face_t* faceFromHe(int32_t heid) { return &faces[halfedges[heid].fid]; }
	const he_t* heFromFace(int32_t fid) const { return &halfedges[faces[fid].heid]; }
	const he_t* heFromVert(int32_t vid) const { return &halfedges[verts[vid].heid]; }
	const vert_t* vertFromHe(int32_t heid) const { return &verts[halfedges[heid].vid]; }
	const face_t* faceFromHe(int32_t heid) const { return &faces[halfedges[heid].fid]; }

	vector<vert_t> verts;
	vector<he_t>   halfedges;
	vector<face_t> faces;
};

HDS_Mesh* buildHalfEdgeMesh(
	const vector<Point3f> &inVerts,
	const vector<PolyIndex> &inFaces
);
void fillNullFaces(
	vector<he_t> &hes,
	vector<face_t> &faces,
	unordered_set<int32_t> &exposedHEs
);