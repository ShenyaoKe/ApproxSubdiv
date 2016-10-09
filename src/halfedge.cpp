#include "halfedge.h"

int32_t HDS_Vertex::uid = 0;
int32_t HDS_HalfEdge::uid = 0;
int32_t HDS_Face::uid = 0;

HDS_Mesh* buildHalfEdgeMesh(
	const vector<Point3f> &inVerts, const vector<PolyIndex> &inFaces)
{
	HDS_Mesh::resetIndex();
	size_t vertsCount = inVerts.size();
	size_t facesCount = inFaces.size();

	size_t heCount = 0;
	// Accumulate face size to get the number of half-edges
	for (size_t i = 0; i < inFaces.size(); i++)
		heCount += inFaces[i].size;

	// Half-Edge arrays for actual HDS_Mesh
	vector<vert_t> verts(vertsCount);
	vector<face_t> faces(facesCount);
	vector<he_t> hes(heCount);
	// Temporary Half-Edge Pair Recorder
	//using hepair_t = pair<int32_t, int32_t>;
	using hepair_t = int64_t;
	auto make_hePair = [](int32_t id1, int32_t id2) {
		return (int64_t(id1) << 32) + int64_t(id2);
	};
	auto reverse_hePair = [](int64_t id) {
		return (id << 32) + (id >> 32);
	};
	// TODO: replace by unordered_map
	unordered_map<hepair_t, int32_t/*, Utils::pair_hash*/> heMap;

	// Assign vertex positions and ids
	for (size_t i = 0; i < vertsCount; i++)
	{
		verts[i].index = i;
	}
	// Malloc Faces
	for (size_t i = 0, heOffset = 0; i < facesCount; i++)
	{
		// Go through all faces
		auto Fi = &inFaces[i];
		int32_t fsize = Fi->size;
		face_t* curFace = &faces[i];

		for (size_t j = 0; j < fsize; j++)
		{
			// calculate current, prev and next edge id
			int32_t curIdx = j + heOffset;

			// link current face and vertex of the edge
			auto &curHe = hes[curIdx];
			// vid in polyindex has offset 1
			curHe.vid = Fi->v[j] - 1;
			curHe.fid = i;
			auto &curVert = verts[curHe.vid];

			// Check index boundary
			// first: prev=last,   next=1
			// last : prev=last-1, next=0
			int32_t jprev = (j == 0) ? fsize - 1 : j - 1;
			int32_t jnext = (j == fsize - 1) ? 0 : j + 1;
			// Connect current edge with previous and next
			curHe.next_offset = jnext - j;
			curHe.prev_offset = jprev - j;

			// connect current vertex to he
			if (curVert.heid == sInvalidHDS) curVert.heid = curHe.index;

			// record edge for flip connection
			int32_t vj = Fi->v[j];
			int32_t vj_next = Fi->v[jnext];
			hepair_t vPair = make_hePair(vj, vj_next);
			// Record edge pair
			if (heMap.find(vPair) == heMap.end())
			{
				heMap[vPair] = curIdx;
			}
		}

		curFace->heid = heOffset;
		//curFace->computeNormal();

		heOffset += Fi->size;
	}
	// hash table for visited edges
	vector<bool> visitedHEs(heMap.size(), false);
	// hash set to record exposed edges
	unordered_set<int32_t> exposedHEs;
	// for each half edge, find its flip
	for (auto heit : heMap)
	{
		int from, to;

		hepair_t hePair = heit.first;
		int32_t heID = heit.second;

		if (!visitedHEs[heID])
		{
			visitedHEs[heID] = true;

			auto invItem = heMap.find(reverse_hePair(hePair));

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

	HDS_Mesh* thismesh = new HDS_Mesh(verts, hes, faces);

	return thismesh;
}

// Functionality: 
//	Add null face and edges directly into original buffer to make mesh validate.
// Input buffers: 
//	half-edges, faces,
//	hash set of indices of exposed edges(flip == null)
// 
void fillNullFaces(
	vector<he_t> &hes,
	vector<face_t> &faces,
	unordered_set<int32_t> &exposedHEs)
{
	// record initial edge number
	size_t initSize = hes.size();
	// allocate unassigned edges
	hes.resize(initSize + exposedHEs.size());
	while (!exposedHEs.empty())
	{
		he_t* he = &hes[*exposedHEs.begin()];

		// Skip checked edges, won't skip in first check
		if (he->flip_offset) continue;

		faces.emplace_back();
		face_t* nullface = &faces.back();

		vector<int32_t> null_hes, null_hefs;
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
			he_t* curHEF = &hes[initSize + i];
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