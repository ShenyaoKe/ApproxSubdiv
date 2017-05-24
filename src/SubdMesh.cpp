#include "SubdMesh.h"
#include "GregroryUtils.h"

SubdMesh::SubdMesh(const char* filename)
{
	if (ObjParser::parse(filename, verts, uvs, norms, fids))
	{
		printf("file %s load successfully!\n", filename);
		mHDSMesh.reset(HDS::buildHalfEdgeMesh(verts, fids));
		initPatch();
		//savePatch();
		evalGregory();
	}
}

SubdMesh::~SubdMesh()
{
}

void SubdMesh::exportIndexedVBO(vector<Float>* vtx_array,
								vector<Float>* /*uv_array*/,
								vector<Float>* /*norm_array*/,
								vector<uint32_t>* idx_array) const
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

void SubdMesh::getPatch(BufferTrait &bezier_trait,
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
	fprintf(fp, "o [BezierPatch]\n");
	uint32_t bezier_patch_size = bezier_patch.size();
	for (int i = 0; i < bezier_patch.size(); i++)
	{
		fprintf(fp, "v %f %f %f\n", bezier_patch[i].x, bezier_patch[i].y, bezier_patch[i].z);
	}
	for (int i = 0; i < bezier_patch.size()/16; i++)
	{
		fprintf(fp, "g [BezierPatchID%04d]\n", i);
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

	fprintf(fp, "o [GregoryPatch]\n");

	for (auto &p : gregory_patch)
	{
		fprintf(fp, "v %f %f %f\n", p.x, p.y, p.z);
	}
	for (int i = 0; i < gregory_patch.size() / sGregoryPatchSize; i++)
	{
		fprintf(fp, "g GregoryPatchID%04d]\n", i);
		int idx = i * sGregoryPatchSize + 1 + bezier_patch_size;
		fprintf(fp,
				// corner
				"f %d %d %d %d %d\n"
				"f %d %d %d %d %d\n"
				"f %d %d %d %d %d\n"
				"f %d %d %d %d %d\n"
				// edge fill
				"f %d %d %d %d\n"
				"f %d %d %d %d\n"
				"f %d %d %d %d\n"
				"f %d %d %d %d\n"
				// face fill
				"f %d %d %d %d %d %d %d %d\n",
				idx + 0,  idx + 1,  idx + 3,  idx + 4,  idx + 2,
				idx + 5,  idx + 6,  idx + 8,  idx + 9,  idx + 7,
				idx + 10, idx + 11, idx + 13, idx + 14, idx + 12,
				idx + 15, idx + 16, idx + 18, idx + 19, idx + 17,

				idx + 1,  idx + 7,  idx + 9,  idx + 3,
				idx + 6,  idx + 12, idx + 14, idx + 8,
				idx + 11, idx + 17, idx + 19, idx + 13,
				idx + 16, idx + 2,  idx + 4,  idx + 18,

				idx + 3,  idx + 9,  idx + 8,  idx + 14,
				idx + 13, idx + 19, idx + 18, idx + 4);
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
	vector<uint32_t> boundaryVertex;

	// Compute vertex valence
	for (int i = 0; i < fids.size() ; i++)
	{
		if (mHDSMesh->faces[i].isNullFace)
		{
			continue;
		}

		auto &vids = fids[i].v;
		bad_face_hash[i] = fids[i].size != 4;

		for_each(vids.begin(), vids.end(), [&vValenceCount] (uint32_t vid) {
			vValenceCount[vid - 1]++;
		});
	}

	// Find all faces using Gregory Patch
	for (int i = 0; i < vValenceCount.size(); i++)
	{
		if (vValenceCount[i] != 4)
		{
			bad_vert_hash[i] = true;
			HDS::HalfEdge* he = mHDSMesh->heFromVert(i);
			HDS::HalfEdge* curHE = he;
			do
			{
				if (curHE->isBoundary)
				{
					boundaryVertex.push_back(i);
				}

				bad_face_hash[curHE->fid] = true;
				curHE = curHE->flip()->next();
			} while (curHE != he);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	genBezierPatch(bad_face_hash);

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
		// Offset in patch array
		if (bad_face_hash[i])
		{
			continue;
		}

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
		HDS::HalfEdge* he = mHDSMesh->heFromFace(i);
		HDS::HalfEdge* curHE = he;
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

		HDS::HalfEdge* he = mHDSMesh->heFromFace(fid);
		HDS::HalfEdge* curHE = he;
		// Loop over current face
		// for each corner
		for (uint32_t j = 0; j < 4; j++)
		{
			uint32_t vid = curHE->vid;
			uint32_t vValence = vValenceCount[vid];
			auto corner_coef = Gregory::corner_coef(vValence);
			Point3f* curP0 = &gregory_patch[pOffset + j * subPatchSize];
			Point3f* curE0_plus = curP0 + 1;
			Point3f* curE0_minus = curP0 + 2;
			// Loop over adjacent edges and vertices around current vertex
			/************************************************************************/
			/* P and e+                                                             */
			/************************************************************************/
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
					*curE0_plus +=
						Gregory::edgeP_coefMi(valenceI, vValence) * heCenters.at(curHEid)
						+ Gregory::edgeP_coefCi(valenceI, vValence) * fCenters.at(curFid);
					
					valenceI++;
				}
				curHE = curHE->rotCCW();
			} while (curHE != he);
			*curP0 *= corner_coef[1];
			*curP0 += corner_coef[0] * verts[vid];
			*curE0_plus *= 4.0 / 3.0 / vValence * Gregory::edge_eigen_val(vValence);
			*curE0_plus += *curP0;
			/************************************************************************/
			/* e-                                                                   */
			/************************************************************************/
			// TODO: remove redundant calculation for e- and p
			valenceI = 0;
			HDS::HalfEdge* prevHE = curHE->rotCCW();
			do
			{
				if (!mHDSMesh->faceFromHe(prevHE->index)->isNullFace)
				{
					// TODO: cache out face center and edge center
					uint32_t curFid = prevHE->fid;
					uint32_t curHEid = prevHE->index;
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
							prevHE->flip()->index, hec
						));
					}
					*curE0_minus += Gregory::edgeP_coefMi(valenceI, vValence) * heCenters.at(curHEid)
						+ Gregory::edgeP_coefCi(valenceI, vValence) * fCenters.at(curFid);

					valenceI++;
				}
				prevHE = prevHE->rotCCW();
			} while (prevHE != curHE->rotCCW());
			*curE0_minus *= 2.0 / vValence;
			*curE0_minus *= 2.0 / 3.0 *  Gregory::edge_eigen_val(vValence);
			*curE0_minus += *curP0;
			// End of e0-

			he = curHE = curHE->next();
		}
		// After all corner points Pi and edge points e+/e- are calculated
		// Calculate face points
		curHE = he;
		for (uint32_t j = 0; j < 4; j++)
		{
			// calculate c0,c1,c2
			uint32_t vid = curHE->vid;
			const int dFaceType = 3;//it's 
			const Float inv_d = 1.0 / dFaceType;
			auto c0 = Gregory::cosPi(2, vValenceCount[vid]);
			auto c1 = Gregory::cosPi(2, vValenceCount[curHE->next()->vid]);
			auto c2 = Gregory::cosPi(2, vValenceCount[curHE->prev()->vid]);

			Point3f* curP0 = &gregory_patch[pOffset + j * subPatchSize];
			Point3f* e_plus = curP0 + 1;
			Point3f* e_minus = curP0 + 2;
			Point3f* f_plus = curP0 + 3;
			Point3f* f_minus = curP0 + 4;
			Point3f* e_prev_plus = j == 0 ? curP0 + 16 : curP0 - 4;
			Point3f* e_next_minus = j == 3 ? curP0 - 13 : curP0 + 7;
			// f+
			// get r+
			Vector3f r_plus =
				(heCenters.at(curHE->prev()->index) - heCenters.at(curHE->rotCW()->index)) / 3.0
				+ (fCenters.at(curHE->fid) - fCenters.at(curHE->rotCW()->fid)) * 2.0 / 3.0;

			*f_plus = (c1 * *curP0
				+ (dFaceType - 2 * c0 - c1) * *e_plus
				+ 2 * c0 * *e_next_minus
				+ r_plus) * inv_d;
			// f-
			// get r-
			Vector3f r_minus =
				(heCenters.at(curHE->index) - heCenters.at(curHE->rotCCW()->prev()->index)) / 3.0
				+ (fCenters.at(curHE->fid) - fCenters.at(curHE->rotCCW()->fid)) * 2.0 / 3.0;

			*f_minus = (c2 * *curP0
				+ (dFaceType - 2 * c0 - c2) * *e_minus
				+ 2 * c0 * *e_prev_plus
				+ r_minus) * inv_d;

			// move to next edge
			curHE = curHE->next();
		}
		pOffset += sGregoryPatchSize;
	}
	gregory_patch.shrink_to_fit();
}

void SubdMesh::evalGregory() const
{
	const int tessSize = 5;
	float us[tessSize];// = { 0, 0.5f, 1.0f };
	float vs[tessSize];// = { 0, 0.5f, 1.0f };
	float invTessSize = 1.0f / tessSize;
	for (int i = 1; i < tessSize - 1; i++)
	{
		us[i] = vs[i] = float(i) * invTessSize;
	}
	us[0] = vs[0] = 0.0f;
	us[tessSize - 1] = vs[tessSize - 1] = 1.0f;

	auto bilinear_pos = [](const Point3f &p0, const Point3f &p1,
						   const Point3f &p2, const Point3f &p3,
						   float s, float t)
	{
		// id3 -- id2
		//  |      |
		// id0 -- id1
		return lerp(lerp(p0, p1, s), lerp(p3, p2, s), t);
	};

	vector<Point3f> ret;
	for (int ui = 0; ui < tessSize; ui++)
	{
		float u = us[ui];
		for (int vi = 0; vi < tessSize; vi++)
		{
			float v = vs[vi];

			if ((u == 0 || u == 1) && (v == 0 || v == 1))
			{
				ret.push_back(bilinear_pos(
					gregory_patch[0],
					gregory_patch[5],
					gregory_patch[10],
					gregory_patch[15],
					u, v
				));
			}
			else
			{
				// Sample on grid
				// Face Points
				Point3f f0 = (u * gregory_patch[3]
							  + v * gregory_patch[4])
					/ (u + v);
				Point3f f1 = ((1 - u) * gregory_patch[9]
							  + v  * gregory_patch[8])
					/ (1 - u + v);
				Point3f f2 = ((1 - u) * gregory_patch[13]
							  + (1 - v) * gregory_patch[14])
					/ (2 - u - v);
				Point3f f3 = (u * gregory_patch[19]
							  + (1 - v) * gregory_patch[18])
					/ (1 + u - v);

				Point3f p[9];
				// 6 -- 7 -- 8
				// |    |    |
				// 3 -- 4 -- 5
				// |    |    |
				// 0 -- 1 -- 2
				p[0] = bilinear_pos(gregory_patch[0],
									gregory_patch[1],
									f0,
									gregory_patch[2],
									u, v);
				p[1] = bilinear_pos(gregory_patch[1],
									gregory_patch[7],
									f1,
									f0,
									u, v);
				p[2] = bilinear_pos(gregory_patch[7],
									gregory_patch[5],
									gregory_patch[6],
									f1,
									u, v);

				p[3] = bilinear_pos(gregory_patch[2],
									f0,
									f3,
									gregory_patch[16],
									u, v);
				p[4] = bilinear_pos(f0, f1, f2, f3,
									u, v);
				p[5] = bilinear_pos(f1,
									gregory_patch[6],
									gregory_patch[12],
									f2,
									u, v);

				p[6] = bilinear_pos(gregory_patch[16],
									f3,
									gregory_patch[17],
									gregory_patch[15],
									u, v);
				p[7] = bilinear_pos(f3,
									f2,
									gregory_patch[11],
									gregory_patch[17],
									u, v);
				p[8] = bilinear_pos(f2,
									gregory_patch[12],
									gregory_patch[10],
									gregory_patch[11],
									u, v);

				// 2 -- 3
				// |    |
				// 0 -- 1
				p[0] = bilinear_pos(p[0], p[1], p[4], p[3], u, v);
				p[1] = bilinear_pos(p[1], p[2], p[5], p[4], u, v);
				p[2] = bilinear_pos(p[3], p[4], p[7], p[6], u, v);
				p[3] = bilinear_pos(p[4], p[5], p[8], p[7], u, v);

				ret.push_back(bilinear_pos(p[0], p[1], p[3], p[2], u, v));
			}
		}
	}
	FILE* fp = fopen("patch_eval.obj", "w");
	for (auto &p : ret)
	{
		fprintf(fp, "v %f %f %f\n", p.x, p.y, p.z);
	}
	for (int i = 0; i < tessSize - 1; i++)
	{
		for (int j = 0; j < tessSize - 1; j++)
		{
			int x = i * tessSize + 1 + j;
			int y = x + tessSize;
			fprintf(fp,
					"f %d %d %d %d\n",
					x, y, y + 1, x + 1);
		}
	}

	fclose(fp);
}
