#include "SubdMesh.h"

SubdMesh::SubdMesh(const char* filename)
{
	if (ObjParser::parse(filename, verts, uvs, norms, fids))
	{
		printf("file %s load successfully!\n", filename);
		mHDSMesh.reset(buildHalfEdgeMesh(verts, fids));
		initPatch();
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
	trait.data = (void*)(patch4.data());
	trait.count = patch4.size();
	trait.size = sizeof(Point3f) * trait.count;
	trait.offset = 0;
	trait.stride = sizeof(Point3f);
}

void SubdMesh::savePatch() const
{
	FILE* fp = fopen("patch.obj", "w");
	for (int i = 0; i < patch4.size(); i++)
	{
		fprintf(fp, "v %f %f %f\n", patch4[i].x, patch4[i].y, patch4[i].z);
	}
	for (int i = 0; i < patch4.size()/16; i++)
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
	fclose(fp);
}

void SubdMesh::initPatch()
{
	vector<uint32_t> vValenceCount(verts.size(), 0);
	vector<bool> bad_faces(fids.size(), false);
	vector<bool> bad_vertices(verts.size(), false);

	for (int i = 0; i < fids.size() ; i++)
	{
		auto &vids = fids[i].v;
		bad_faces[i] = fids[i].size != 4;

		for_each(vids.begin(), vids.end(), [&vValenceCount] (uint32_t vid) {
			vValenceCount[vid - 1]++;
		});
	}
	uint32_t badvcount = 0;
	for (int i = 0; i < vValenceCount.size(); i++)
	{
		if (vValenceCount[i] != 4)
		{
			bad_vertices[i] = true;
			badvcount++;
			auto he = mHDSMesh->heFromVert(i);
			auto curHE = he;
			do 
			{
				bad_faces[curHE->fid] = true;
				curHE = curHE->flip()->next();
			} while (curHE != he);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	// Patch Generation
	patch4.clear();

	patch4.reserve(fids.size() * sPatchSize);
	for (size_t i = 0, pOffset = 0; i < fids.size(); i++)
	{
		auto &fid = fids[i];
		// Offset in patch array
		if (bad_faces[i]) continue;

		patch4.resize(pOffset + sPatchSize);
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
		patch4[pOffset] = cornerPatch(5);
		patch4[pOffset + 3] = cornerPatch(6);
		patch4[pOffset + 12] = cornerPatch(9);
		patch4[pOffset + 15] = cornerPatch(10);

		patch4[pOffset + 1] = edgePatch(5, 6);
		patch4[pOffset + 2] = edgePatch(6, 5);
		patch4[pOffset + 4] = edgePatch(5, 9);
		patch4[pOffset + 8] = edgePatch(9, 5);
		patch4[pOffset + 7] = edgePatch(6, 10);
		patch4[pOffset + 11] = edgePatch(10, 6);
		patch4[pOffset + 13] = edgePatch(9, 10);
		patch4[pOffset + 14] = edgePatch(10, 9);

		patch4[pOffset + 5] = interialPatch(5, 6, 9, 10);
		patch4[pOffset + 6] = interialPatch(6, 5, 10, 9);
		patch4[pOffset + 9] = interialPatch(9, 5, 10, 6);
		patch4[pOffset + 10] = interialPatch(10, 6, 9, 5);
		// Move on
		pOffset += sPatchSize;
	}
}
