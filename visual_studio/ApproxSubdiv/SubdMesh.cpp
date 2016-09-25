#include "SubdMesh.h"

SubdMesh::SubdMesh(const char* filename)
{
	if (ObjParser::parse(filename, verts, uvs, norms, fids))
	{
		printf("file %s load successfully!\n", filename);
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
