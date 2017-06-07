#include "HalfEdge.h"

namespace HDS
{

SizeType Vertex::uid = 0;
SizeType HalfEdge::uid = 0;
SizeType Face::uid = 0;

Mesh* Mesh::buildMesh(SizeType vertexCount,
					  const vector<SizeType> &faceIds,
					  const vector<SizeType> &faceSide,
					  const vector<SizeType> &faceIdOffset)
{
	Mesh::resetIndex();

	constexpr OffsetType idOfs = 0;

	size_t facesCount = faceIdOffset.size();

	size_t heCount = 0;
	// Accumulate face size to get the number of half-edges
	for (size_t i = 0; i < facesCount; i++)
	{
		heCount += faceSide[i];
	}

	// Half-Edge arrays for actual Mesh
	vector<Vertex> verts(vertexCount);
	vector<Face> faces(facesCount);
	vector<HalfEdge> hes(heCount);
	// Temporary Half-Edge Pair Recorder
	//using hepair_t = pair<SizeType, SizeType>;
	using hepair_t = int64_t;
	auto make_hePair = [] (SizeType id1, SizeType id2) {
		return (int64_t(id1) << 32) + int64_t(id2);
	};
	auto reverse_hePair = [] (int64_t id) {
		return (id << 32) + (id >> 32);
	};

	unordered_map<hepair_t, SizeType> heMap;

	// Malloc Faces
	for (size_t i = 0, heOffset = 0; i < facesCount; i++)
	{
		// Go through all faces
		const SizeType* fids = &faceIds[faceIdOffset[i]];
		SizeType curFaceSide = faceSide[i];
		Face* curFace = &faces[i];

		for (size_t j = 0; j < curFaceSide; j++)
		{
			// calculate current, prev and next edge id
			SizeType curIdx = j + heOffset;

			// link current face and vertex of the edge
			auto &curHe = hes[curIdx];
			// vid from obj index has offset 1
			curHe.vid = fids[j] + idOfs;
			curHe.fid = i;
			auto &curVert = verts[curHe.vid];

			// Check index boundary
			// first: prev=last,   next=1
			// last : prev=last-1, next=0
			SizeType jprev = (j == 0) ? curFaceSide - 1 : j - 1;
			SizeType jnext = (j == curFaceSide - 1) ? 0 : j + 1;
			// Connect current edge with previous and next
			curHe.next_offset = jnext - j;
			curHe.prev_offset = jprev - j;

			// connect current vertex to he
			if (curVert.heid == cInvalidIndex)
			{
				curVert.heid = curHe.index;
			}

			// record edge for flip connection
			SizeType vj = fids[j];
			SizeType vj_next = fids[jnext];
			hepair_t vPair = make_hePair(vj, vj_next);
			// Record edge pair
			if (heMap.find(vPair) == heMap.end())
			{
				heMap[vPair] = curIdx;
			}
		}

		curFace->heid = heOffset;

		heOffset += curFaceSide;
	}
	// hash table for visited edges
	vector<bool> visitedHEs(heCount, false);
	// hash set to record exposed edges
	unordered_set<SizeType> exposedHEs;
	// for each half edge, find its flip
	for (auto &heit : heMap)
	{
		hepair_t hePair = heit.first;
		SizeType heID = heit.second;

		if (!visitedHEs[heID])
		{
			visitedHEs[heID] = true;

			auto &invItem = heMap.find(reverse_hePair(hePair));

			if (invItem != heMap.end())
			{
				auto &he = hes[heit.second];
				auto &hef = hes[invItem->second];

				he.flip_offset = hef.index - he.index;
				hef.flip_offset = -he.flip_offset;

				visitedHEs[invItem->second] = true;
			}
			else
			{
				exposedHEs.insert(heID);
			}
		}
	}

	// Check Holes and Fill with Null Faces
	if (!exposedHEs.empty())
	{
		// fill exposed edges (flip == null) with null faces
		fillNullFaces(hes, faces, exposedHEs);
	}

	Mesh* thismesh = new Mesh(verts, hes, faces);

	return thismesh;
}

void fillNullFaces(vector<HalfEdge> &hes,
				   vector<Face> &faces,
				   unordered_set<SizeType> &exposedHEs)
{
	// record initial edge number
	size_t initSize = hes.size();
	// allocate unassigned edges
	hes.resize(initSize + exposedHEs.size());
	while (!exposedHEs.empty())
	{
		HalfEdge* he = &hes[*exposedHEs.begin()];

		// Skip checked edges, won't skip in first check
		if (he->flip_offset) continue;

		faces.emplace_back();
		Face* nullface = &faces.back();

		vector<SizeType> null_hes, null_hefs;
		size_t heIdOffset = initSize;
		// record null edges on the same null face
		auto curHE = he;
		do
		{
			null_hes.push_back(curHE->index);
			exposedHEs.erase(curHE->index);
			null_hefs.push_back(heIdOffset++);

			// if curHE->next->flip == null (offset != 0),
			//     found the next exposed edge
			///                       ___curHE___
			///                      /
			///     exposed edge--> / curHE->next
			///                    /
			// else, move to curHE->next->flip->next
			///                    \    <--exposed edge
			///         curHE->nex  \ 
			///         ->flip->next \  ___curHE___ 
			///                      / / 
			///         curHE->next / /curHE->next
			///         ->flip     / / 
			curHE = curHE->next();
			// Loop adjacent edges to find the exposed edge
			while (curHE->flip_offset)
			{
				curHE = curHE->flip()->next();
			}
		} while (curHE != he);
		// get edge number of current null face
		size_t nNullEdges = null_hes.size();

		// construct null face
		for (size_t i = 0; i < nNullEdges; i++)
		{
			curHE = &hes[null_hes[i]];
			HalfEdge* curHEF = &hes[initSize + i];
			null_hefs[i] = initSize + i;

			curHE->isBoundary = curHEF->isBoundary = true;
			curHE->flip_offset = curHEF - curHE;
			curHEF->flip_offset = -curHE->flip_offset;
			curHEF->vid = curHE->next()->vid;
			curHEF->fid = nullface->index;

			/// Buffer: ...(existing edges)..., 0, 1, 2, ..., n-1
			/// Structure:   e(n-1)-> ... -> e1 -> e0 -> e(n-1)
			// prev edge is the next one in buffer,
			// except the last one, previous edge is the first one in buffer
			curHEF->prev_offset = (i < nNullEdges - 1) ? 1 : 1 - nNullEdges;
			// next edge is the previous one in buffer,
			// except the first one, next edge is the last one in buffer
			curHEF->next_offset = (i > 0) ? -1 : nNullEdges - 1;
		}
		// Update Null Face Component and Flag
		nullface->isNullFace = true;
		//nullface->heid = he->index;
		nullface->heid = null_hefs[0];

		initSize += nNullEdges;
	}
}

}
