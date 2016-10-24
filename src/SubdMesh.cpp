#include "SubdMesh.h"
#include "GregroryPatch.hh"

SubdMesh::SubdMesh(const char* filename)
{
	if (ObjParser::parse(filename, verts, uvs, norms, fids))
	{
		printf("file %s load successfully!\n", filename);
		mHDSMesh.reset(buildHalfEdgeMesh(verts, fids));
		initPatch();
		savePatch();
	}
}

SubdMesh::~SubdMesh()
{
}

void SubdMesh::exportIndexedVBO(
	vector<Float>* vtx_array,
	vector<Float>* uv_array,
	vector<Float>* norm_array,
	vector<unsigned int>* idx_array) const
{
	vtx_array->reserve(verts.size() * 3);
	for (int i = 0; i < verts.size(); i++)
	{
		vtx_array->push_back(verts[i].x);
		vtx_array->push_back(verts[i].y);
		vtx_array->push_back(verts[i].z);
	}
	idx_array->reserve(fids.size() * 2);
	for (int i = 0; i < fids.size(); i++)
	{
		auto &fid = fids[i];
		if (fid.size == 4)
		{
			idx_array->push_back(fid.v[0] - 1);
			idx_array->push_back(fid.v[1] - 1);
			idx_array->push_back(fid.v[2] - 1);
			idx_array->push_back(fid.v[3] - 1);
		}
	}
}

void SubdMesh::getPatch(BufferTrait &trait) const
{
	trait.data = (void*)(bezier_patch.data());
	trait.count = bezier_patch.size();
	trait.size = sizeof(Point3f) * trait.count;
	trait.offset = 0;
	trait.stride = sizeof(Point3f);
}

void SubdMesh::getPatch(
	BufferTrait &bezier_trait,
	BufferTrait &gregory_trait) const
{
	bezier_trait.data = (void*)(bezier_patch.data());
	bezier_trait.count = bezier_patch.size();
	bezier_trait.size = sizeof(Point3f) * bezier_trait.count;
	bezier_trait.offset = 0;
	bezier_trait.stride = sizeof(Point3f);

	gregory_trait.data = (void*)(gregory_patch.data());
	gregory_trait.count = gregory_patch.size();
	gregory_trait.size = sizeof(Point3f) * gregory_trait.count;
	gregory_trait.offset = 0;
	gregory_trait.stride = sizeof(Point3f);
}

void SubdMesh::savePatch() const
{
	FILE* fp = fopen("patch.obj", "w");
	uint32_t bezier_patch_size = bezier_patch.size();
	for (int i = 0; i < bezier_patch.size(); i++)
	{
		fprintf(fp, "v %f %f %f\n", bezier_patch[i].x, bezier_patch[i].y, bezier_patch[i].z);
	}
	for (int i = 0; i < bezier_patch.size()/16; i++)
	{
		int idx = i * 16 + 1;
		fprintf(fp, "f %d %d %d %d\n", idx + 0, idx + 1, idx + 5, idx + 4);
		fprintf(fp, "f %d %d %d %d\n", idx + 1, idx + 2, idx + 6, idx + 5);
		fprintf(fp, "f %d %d %d %d\n", idx + 2, idx + 3, idx + 7, idx + 6);

		fprintf(fp, "f %d %d %d %d\n", idx + 4, idx + 5, idx + 9, idx + 8);
		fprintf(fp, "f %d %d %d %d\n", idx + 5, idx + 6, idx + 10, idx + 9);
		fprintf(fp, "f %d %d %d %d\n", idx + 6, idx + 7, idx + 11, idx + 10);

		fprintf(fp, "f %d %d %d %d\n", idx + 8, idx + 9, idx + 13, idx + 12);
		fprintf(fp, "f %d %d %d %d\n", idx + 9, idx + 10, idx + 14, idx + 13);
		fprintf(fp, "f %d %d %d %d\n", idx + 10, idx + 11, idx + 15, idx + 14);
	}
	for (auto &p : gregory_patch)
	{
		fprintf(fp, "v %f %f %f\n", p.x, p.y, p.z);
	}
	for (int i = 0; i < gregory_patch.size() / sGregoryPatchSize; i++)
	{
		int idx = i * sGregoryPatchSize + 1 + bezier_patch_size;
		fprintf(fp, "f %d %d %d %d %d %d %d %d\n",
			idx + 0, idx + 1,
			idx + 5, idx + 6,
			idx + 10, idx + 11,
			idx + 15, idx + 16);
	}
	fclose(fp);
}

Point3f SubdMesh::faceCenter(uint32_t fid) const
{
	auto he = mHDSMesh->heFromFace(fid);
	auto curHE = he;
	Point3f cp;
	uint32_t pCount = 0;
	do
	{
		cp += verts[curHE->vid];
		pCount++;
		curHE = curHE->next();
	} while (curHE != he);
	if (pCount > 0)
	{
		cp /= pCount;
	}
	return cp;
}

Point3f SubdMesh::edgeCenter(uint32_t heid) const
{
	auto he = &mHDSMesh->halfedges[heid];
	return (verts[he->vid] + verts[he->flip()->vid]) * 0.5f;
}

void SubdMesh::initPatch()
{
	vector<uint32_t> vValenceCount(verts.size(), 0);
	vector<uint32_t> badfaces;
	vector<bool> bad_face_hash(fids.size(), false);
	vector<bool> bad_vert_hash(verts.size(), false);

	for (int i = 0; i < fids.size() ; i++)
	{
		auto &vids = fids[i].v;
		bad_face_hash[i] = fids[i].size != 4;

		for_each(vids.begin(), vids.end(), [&vValenceCount] (uint32_t vid) {
			vValenceCount[vid - 1]++;
		});
	}
	uint32_t badvcount = 0;
	for (int i = 0; i < vValenceCount.size(); i++)
	{
		if (vValenceCount[i] != 4)
		{
			bad_vert_hash[i] = true;
			badvcount++;
			he_t* he = mHDSMesh->heFromVert(i);
			he_t* curHE = he;
			do
			{
				bad_face_hash[curHE->fid] = true;
				curHE = curHE->flip()->next();
			} while (curHE != he);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	genBezierPatch(bad_face_hash);

	int badcount = 0;
	for (int i = 0; i < bad_face_hash.size(); i++)
	{
		if (bad_face_hash[i])
		{
			badfaces.push_back(i);
		}
	}
	genGregoryPatch(vValenceCount, badfaces);
}

void SubdMesh::genBezierPatch(const vector<bool> &bad_face_hash)
{
	// Patch Generation
	bezier_patch.clear();

	bezier_patch.reserve(fids.size() * sBezierPatchSize);
	for (size_t i = 0, pOffset = 0; i < fids.size(); i++)
	{
		auto &fid = fids[i];
		// Offset in patch array
		if (bad_face_hash[i]) continue;

		bezier_patch.resize(pOffset + sBezierPatchSize);
		// calculate patch
		uint32_t vid[16];
		//////////////////////////////////////////////////////////////////////////
		// Cache out neighbour vertices
		auto he = mHDSMesh->heFromFace(i);
		auto he_oppo = he->flip()->prev();
		vid[5] = he->vid;
		vid[1] = he_oppo->prev()->vid;
		vid[2] = he_oppo->vid;
		vid[3] = he_oppo->flip()->next()->next()->vid;

		he = he->next();
		he_oppo = he->flip()->prev();
		vid[6] = he->vid;
		vid[7] = he_oppo->prev()->vid;
		vid[11] = he_oppo->vid;
		vid[15] = he_oppo->flip()->next()->next()->vid;

		he = he->next();
		he_oppo = he->flip()->prev();
		vid[10] = he->vid;
		vid[14] = he_oppo->prev()->vid;
		vid[13] = he_oppo->vid;
		vid[12] = he_oppo->flip()->next()->next()->vid;

		he = he->next();
		he_oppo = he->flip()->prev();
		vid[9] = he->vid;
		vid[8] = he_oppo->prev()->vid;
		vid[4] = he_oppo->vid;
		vid[0] = he_oppo->flip()->next()->next()->vid;
		//////////////////////////////////////////////////////////////////////////

		auto cornerPatch = [&](uint32_t id) {
			return (verts[vid[id]] * 16
				+ (verts[vid[id - 1]] + verts[vid[id + 1]]
					+ verts[vid[id - 4]] + verts[vid[id + 4]]) * 4
				+ verts[vid[id - 5]] + verts[vid[id + 5]]
				+ verts[vid[id - 3]] + verts[vid[id + 3]]) / 36.0f;
		};
		auto edgePatch = [&](uint32_t id, uint32_t id2) {
			auto neighbor_off = 5u - abs(int32_t(id2 - id));
			return (
				verts[vid[id]] * 8
				+ verts[vid[id2]] * 4
				+ (verts[vid[id + neighbor_off]] + verts[vid[id - neighbor_off]]) * 2
				+ verts[vid[id2 + neighbor_off]] + verts[vid[id2 - neighbor_off]]
				) / 18.0f;
		};
		auto interialPatch = [&]
		(uint32_t id, uint32_t id2, uint32_t id3, uint32_t id_diag) {
			return (verts[vid[id]] * 4
				+ (verts[vid[id2]] + verts[vid[id3]]) * 2
				+ verts[vid[id_diag]]) / 9.0f;
		};
		bezier_patch[pOffset] = cornerPatch(5);
		bezier_patch[pOffset + 3] = cornerPatch(6);
		bezier_patch[pOffset + 12] = cornerPatch(9);
		bezier_patch[pOffset + 15] = cornerPatch(10);

		bezier_patch[pOffset + 1] = edgePatch(5, 6);
		bezier_patch[pOffset + 2] = edgePatch(6, 5);
		bezier_patch[pOffset + 4] = edgePatch(5, 9);
		bezier_patch[pOffset + 8] = edgePatch(9, 5);
		bezier_patch[pOffset + 7] = edgePatch(6, 10);
		bezier_patch[pOffset + 11] = edgePatch(10, 6);
		bezier_patch[pOffset + 13] = edgePatch(9, 10);
		bezier_patch[pOffset + 14] = edgePatch(10, 9);

		bezier_patch[pOffset + 5] = interialPatch(5, 6, 9, 10);
		bezier_patch[pOffset + 6] = interialPatch(6, 5, 10, 9);
		bezier_patch[pOffset + 9] = interialPatch(9, 5, 10, 6);
		bezier_patch[pOffset + 10] = interialPatch(10, 6, 9, 5);
		// Move on
		pOffset += sBezierPatchSize;
	}
	bezier_patch.shrink_to_fit();
}

void SubdMesh::genGregoryPatch(
	vector<uint32_t> &vValenceCount,
	vector<uint32_t> &irreg_faces
)
{
	unordered_map<uint32_t, Point3f> fCenters;
	unordered_map<uint32_t, Point3f> heCenters;
	uint32_t patchCount = irreg_faces.size();
	// calculate face center and edge center
	for (uint32_t i = 0; i < patchCount; i++)
	{
		he_t* he = mHDSMesh->heFromFace(i);
		he_t* curHE = he;
		Point3f fc(0, 0, 0);
		do
		{
			if (heCenters.find(curHE->index) == heCenters.end())
			{
				auto hec = (verts[curHE->vid]
					+ verts[curHE->flip()->vid]) * 0.5f;
				heCenters.insert(make_pair(curHE->index, hec));
				heCenters.insert(make_pair(curHE->flip()->index, hec));
			}
			fc += verts[curHE->vid];
			curHE = curHE->next();
		} while (curHE != he);
		fCenters.insert(make_pair(i, fc * 0.25f));
	}

	gregory_patch.reserve(patchCount * sGregoryPatchSize);

	size_t pOffset = 0;
	const uint32_t subPatchSize = 5;
	// Go through faces
	for (uint32_t i = 0; i < patchCount; i++)
	{
		uint32_t fid = irreg_faces[i];
		PolyIndex &fIdx = fids[fid];
		if (fIdx.size != 4) continue;

		gregory_patch.resize(pOffset + sGregoryPatchSize);

		he_t* he = mHDSMesh->heFromFace(fid);
		he_t* curHE = he;
		// Loop over current face
		for (uint32_t j = 0; j < 4; j++)
		{
			uint32_t vid = curHE->vid;
			uint32_t vValence = vValenceCount[vid];
			auto corner_coef = Gregory::corner_coef(vValence);
			Point3f* curP0 = &gregory_patch[pOffset + j * subPatchSize];
			Point3f* curE0_plus = curP0 + 1;
			// Loop over current vertex
			uint32_t valenceI = 0;
			do
			{
				if (!mHDSMesh->faceFromHe(curHE->index)->isNullFace)
				{
					// TODO: cache out face center and edge center
					uint32_t curFid = curHE->fid;
					uint32_t curHEid = curHE->index;
					if (fCenters.find(curFid) == fCenters.end())
					{
						fCenters.insert(make_pair(
							curFid, faceCenter(curFid)
						));
					}
					if (heCenters.find(curHEid) == heCenters.end())
					{
						auto hec = edgeCenter(curHEid);
						heCenters.insert(make_pair(
							curHEid, hec
						));
						heCenters.insert(make_pair(
							curHE->flip()->index, hec
						));
					}
					// Accumulate corner point
					*curP0 += fCenters.at(curFid) + heCenters.at(curHEid);
					*curE0_plus += Gregory::edgeP_coefMi(valenceI, vValence) * heCenters.at(curHEid)
						+ Gregory::edgeP_coefCi(valenceI, vValence) * fCenters.at(curFid);

					valenceI++;
				}
				curHE = curHE->rotCCW();
			} while (curHE != he);
			*curP0 *= corner_coef[1];
			*curP0 += corner_coef[0] * verts[vid];
			*curE0_plus *= 4.0 / 3.0 / vValence * Gregory::edge_lambda(vValence);
			*curE0_plus += *curP0;

			he = curHE = curHE->next();
		}
		pOffset += sGregoryPatchSize;
	}
	gregory_patch.shrink_to_fit();
}
