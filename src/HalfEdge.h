#pragma once
#include "common.h"
#include "Geometry/TriangleMesh.h"

using SizeType = uint32_t;
using OffsetType = int32_t;

static const SizeType cInvalidIndex = static_cast<SizeType>(-1);

namespace HDS
{

// Vertex
class Vertex
{
public:
	Vertex() : index(uid++), pid(cInvalidIndex), heid(cInvalidIndex) {}
	~Vertex() {}

	static void resetIndex() { uid = 0; }

public:
	SizeType pid;
	SizeType index;
	SizeType heid;

private:
	static SizeType uid;
};


// Half-Edge
class HalfEdge
{
public:
	static void resetIndex() { uid = 0; }
	static void matchIndexToSize(SizeType size) { uid = size; }

	HalfEdge() : index(uid++), fid(cInvalidIndex), vid(cInvalidIndex)
		, prev_offset(0), next_offset(0), flip_offset(0)
	{
	}
	~HalfEdge() {}

	// Get the explicit pointer to corresponding edges
	HalfEdge* prev() { return this + prev_offset; }
	HalfEdge* next() { return this + next_offset; }
	HalfEdge* flip() { return this + flip_offset; }
	HalfEdge* rotCW() { return flip()->next(); }
	HalfEdge* rotCCW() { return prev()->flip(); }
	const HalfEdge* prev() const { return this + prev_offset; }
	const HalfEdge* next() const { return this + next_offset; }
	const HalfEdge* flip() const { return this + flip_offset; }
	const HalfEdge* rotCW() const { return flip()->next(); }
	const HalfEdge* rotCCW() const { return prev()->flip(); }

	//bool isBoundary() const { return flip_offset == 0; }

	// No self-loop edge
	void setToInvalid() { prev_offset = next_offset = 0; fid = cInvalidIndex; }
	bool isInvalid() const { return prev_offset == 0 || next_offset == 0 || fid == cInvalidIndex; }

	void setFlip(HalfEdge* f_e)
	{
		flip_offset = f_e - this;
		f_e->flip_offset = -flip_offset;
	}
	void breakFlip()
	{
		flip_offset = flip()->flip_offset = 0;
	}

public:
	SizeType index;
	SizeType fid;
	SizeType vid;
	// Offset to index of previous/next/flip edge
	// previous/next/flip edge doesn't exist
	// when (previous/nex/flip == 0)
	OffsetType prev_offset, next_offset, flip_offset;

	bool isBoundary = false;

private:
	static SizeType uid;
};

// Face
class Face
{
public:
	static void resetIndex() { uid = 0; }

	Face() : index(uid++), heid(cInvalidIndex) {}
	~Face() {}

	// Get the connected half-edge id
	// Explicit pointer access is handled by HDS_Mesh
	SizeType heID() const { return heid; }


	void setToInvalid() { heid = cInvalidIndex; }
	bool isInvalid() const { return heid == cInvalidIndex; }

public:
	// Member data
	SizeType index;
	SizeType heid;
	bool isNullFace = false;
private:
	static SizeType uid;
};

// Mesh
class Mesh
{
public:
	Mesh() {}
	Mesh(vector<Vertex> &vs, vector<HalfEdge> &hes, vector<Face> &fs)
	: verts(std::move(vs))
	, halfedges(std::move(hes))
	, faces(std::move(fs)) {}
	Mesh(const Mesh &other)
		: verts(other.verts), halfedges(other.halfedges), faces(other.faces)
	{
	}
	~Mesh() {}

	// Reset UID in each component
	// Mask(Bitwise Operation): face|edge|vertex
	//    e.g. All(111==3), Vertex Only(001==1), Vertex+Edge(011==3)
	static void resetIndex(uint8_t reset_mask = 7)
	{
		if (reset_mask & 1) Vertex::resetIndex();
		if (reset_mask & 2) HalfEdge::resetIndex();
		if (reset_mask & 4) Face::resetIndex();
	}


	void printInfo(const std::string &msg = "")
	{
		if (!msg.empty())
		{
			std::cout << msg << std::endl;
		}
		std::cout << "#vertices = " << verts.size() << std::endl;
		std::cout << "#faces = " << faces.size() << std::endl;
		std::cout << "#half edges = " << halfedges.size() << std::endl;
	}

	HalfEdge* heFromFace(SizeType fid) { return &halfedges[faces[fid].heid]; }
	HalfEdge* heFromVert(SizeType vid) { return &halfedges[verts[vid].heid]; }
	Vertex* vertFromHe(SizeType heid) { return &verts[halfedges[heid].vid]; }
	Face* faceFromHe(SizeType heid) { return &faces[halfedges[heid].fid]; }
	const HalfEdge* heFromFace(SizeType fid) const { return &halfedges[faces[fid].heid]; }
	const HalfEdge* heFromVert(SizeType vid) const { return &halfedges[verts[vid].heid]; }
	const Vertex* vertFromHe(SizeType heid) const { return &verts[halfedges[heid].vid]; }
	const Face* faceFromHe(SizeType heid) const { return &faces[halfedges[heid].fid]; }

public:
	vector<Vertex> verts;
	vector<HalfEdge> halfedges;
	vector<Face>     faces;
};

static Mesh* buildHalfEdgeMesh(SizeType vertexCount,
							   const vector<SizeType> &faceIds,
							   const vector<SizeType> &faceSide,
							   const vector<SizeType> &faceIdOffset);
// Functionality: 
//     Add null face and edges directly into original buffer to make mesh validate.
// Input buffers: 
//     half-edges, faces,
//     hash set of indices of exposed edges(flip == null)
static void fillNullFaces(vector<HalfEdge> &hes,
						  vector<Face> &faces,
						  unordered_set<SizeType> &exposedHEs);
}
