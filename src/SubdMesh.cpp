#include "SubdMesh.h"

SubdMesh::SubdMesh(const char* filename)
{
	if (ObjParser::parse(filename, verts, uvs, norms, fids))
	{
		printf("file %s load successfully!\n", filename);
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
/*

			idx_array->push_back(fid.v[0] - 1);
			idx_array->push_back(fid.v[2] - 1);
			idx_array->push_back(fid.v[3] - 1);*/
		}
		/*else if (fid.size == 3)
		{
			idx_array->push_back(fid.v[0]);
			idx_array->push_back(fid.v[1]);
			idx_array->push_back(fid.v[2]);
		}*/
	}
}

void SubdMesh::getPatch(BufferTrait &trait) const
{
	trait.data = (void*)(&patch[0].x);
	trait.size = sizeof(Point3f) * verts.size();
	trait.offset = 0;
	trait.stride = sizeof(Point3f);
}

void SubdMesh::initPatch()
{
	patch.clear();

	patch.resize(fids.size() * sPatchSize);
	for (size_t i = 0; i < fids.size(); i++)
	{
		auto &fid = fids[i];
		// Offset in patch array
		size_t subId = i * sPatchSize;
		if (fid.size == 4)
		{
			patch[subId     ] = verts[fid.v[0]];
			patch[subId +  3] = verts[fid.v[1]];
			patch[subId + 12] = verts[fid.v[3]];
			patch[subId + 15] = verts[fid.v[2]];

			patch[subId +  1] = lerp(patch[subId], patch[subId + 3], 0.333f);
			patch[subId +  2] = lerp(patch[subId], patch[subId + 3], 0.667f);
			patch[subId +  4] = lerp(patch[subId], patch[subId + 12], 0.333f);
			patch[subId +  8] = lerp(patch[subId], patch[subId + 12], 0.667f);
			patch[subId +  7] = lerp(patch[subId +  3], patch[subId + 15], 0.333f);
			patch[subId + 11] = lerp(patch[subId +  3], patch[subId + 15], 0.667f);
			patch[subId + 13] = lerp(patch[subId + 12], patch[subId + 15], 0.333f);
			patch[subId + 14] = lerp(patch[subId + 12], patch[subId + 15], 0.667f);


			patch[subId +  5] = lerp(patch[subId + 4], patch[subId +  7], 0.333f);
			patch[subId +  6] = lerp(patch[subId + 4], patch[subId +  7], 0.667f);
			patch[subId +  9] = lerp(patch[subId + 8], patch[subId + 11], 0.333f);
			patch[subId + 10] = lerp(patch[subId + 8], patch[subId + 11], 0.667f);
		}
	}
}
