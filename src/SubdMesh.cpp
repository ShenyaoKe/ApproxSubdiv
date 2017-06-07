#include "SubdMesh.h"
#include "GregroryUtils.h"
#include "Kaguya/IO/ObjLoader.h"

namespace Kaguya
{

SubdMesh::SubdMesh(const std::string &filename)
{
	ObjBuffers tmpBuffers;
	if (ObjLoader::loadRawBuffers(tmpBuffers, filename))
	{
		mVerts = std::move(tmpBuffers.vertexBuffer);
		mFaceIndices = std::move(tmpBuffers.faceIndexBuffer);
		mFaceSideCount = std::move(tmpBuffers.faceSizeBuffer);
		mFaceIdOffset.resize(mFaceSideCount.size(), 0);
		for (SizeType i = 0; i < mFaceIdOffset.size() - 1; i++)
		{
			mFaceIdOffset[i + 1] = mFaceIdOffset[i] + mFaceSideCount[i];
		}

		printf("file %s load successfully!\n", filename.c_str());
		mHDSMesh.reset(HDS::Mesh::buildMesh(mVerts.size(),
										mFaceIndices,
										mFaceSideCount,
										mFaceIdOffset));
		process();
#if 0
		savePatch();
#endif
		//evalGregory();
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
	vtx_array->reserve(mVerts.size() * 3);
	for (int i = 0; i < mVerts.size(); i++)
	{
		vtx_array->push_back(mVerts[i].x);
		vtx_array->push_back(mVerts[i].y);
		vtx_array->push_back(mVerts[i].z);
	}
	idx_array->reserve(mFaceIndices.size());
	for (int i = 0; i < mFaceIdOffset.size(); i++)
	{
		const SizeType faceSide = mFaceSideCount[i];
		const SizeType* vids = &mFaceIndices[mFaceIdOffset[i]];
		if (faceSide == 4)
		{
			idx_array->push_back(vids[0] - 1);
			idx_array->push_back(vids[1] - 1);
			idx_array->push_back(vids[2] - 1);
			idx_array->push_back(vids[3] - 1);
		}
	}
}

void SubdMesh::getPatch(BufferTrait &vertexBuffer,
						BufferTrait &bezierPatchIndexBuffer,
						BufferTrait &quadGregoryPatchIndexBuffer,
						BufferTrait &triGregoryPatchIndexBuffer) const
{
	vertexBuffer.data = (void*)mPatchVertexBuffer.data();
	vertexBuffer.count = mPatchVertexBuffer.size();
	vertexBuffer.size = sizeof(Point3f) * vertexBuffer.count;
	vertexBuffer.offset = 0;
	vertexBuffer.stride = sizeof(Point3f);

	bezierPatchIndexBuffer.data = (void*)(mBezierPatchIndices.data());
	bezierPatchIndexBuffer.count = mBezierPatchIndices.size();
	bezierPatchIndexBuffer.size = sizeof(SizeType) * bezierPatchIndexBuffer.count;
	bezierPatchIndexBuffer.offset = 0;
	bezierPatchIndexBuffer.stride = sizeof(SizeType);

	quadGregoryPatchIndexBuffer.data = (void*)(mQuadGregoryPatchIndices.data());
	quadGregoryPatchIndexBuffer.count = mQuadGregoryPatchIndices.size();
	quadGregoryPatchIndexBuffer.size = sizeof(SizeType) * quadGregoryPatchIndexBuffer.count;
	quadGregoryPatchIndexBuffer.offset = 0;
	quadGregoryPatchIndexBuffer.stride = sizeof(SizeType);

	triGregoryPatchIndexBuffer.data = (void*)(mTriGregoryPatchIndices.data());
	triGregoryPatchIndexBuffer.count = mTriGregoryPatchIndices.size();
	triGregoryPatchIndexBuffer.size = sizeof(SizeType) * triGregoryPatchIndexBuffer.count;
	triGregoryPatchIndexBuffer.offset = 0;
	triGregoryPatchIndexBuffer.stride = sizeof(SizeType);
}

void SubdMesh::savePatch() const
{
	FILE* fp = fopen("patch.obj", "w");
	SizeType vertexSize = mPatchVertexBuffer.size();
	for (int i = 0; i < vertexSize; i++)
	{
		fprintf(fp, "v %f %f %f\n",
				mPatchVertexBuffer[i].x,
				mPatchVertexBuffer[i].y,
				mPatchVertexBuffer[i].z);
	}

	fprintf(fp, "o [BezierPatch]\n");
	for (int i = 0; i < mBezierPatchIndices.size() / sQuadBezierPatchSize; i++)
	{
		fprintf(fp, "g [BezierPatchID%04d]\n", i);
		const SizeType* idx = &mBezierPatchIndices[i * 16];
		fprintf(fp, "f %d %d %d %d\n", idx[0] + 1, idx[1] + 1, idx[5] + 1, idx[4] + 1);
		fprintf(fp, "f %d %d %d %d\n", idx[1] + 1, idx[2] + 1, idx[6] + 1, idx[5] + 1);
		fprintf(fp, "f %d %d %d %d\n", idx[2] + 1, idx[3] + 1, idx[7] + 1, idx[6] + 1);

		fprintf(fp, "f %d %d %d %d\n", idx[4] + 1, idx[5] + 1, idx[9] + 1, idx[8] + 1);
		fprintf(fp, "f %d %d %d %d\n", idx[5] + 1, idx[6] + 1, idx[10] + 1, idx[9] + 1);
		fprintf(fp, "f %d %d %d %d\n", idx[6] + 1, idx[7] + 1, idx[11] + 1, idx[10] + 1);

		fprintf(fp, "f %d %d %d %d\n", idx[8] + 1, idx[9] + 1, idx[13] + 1, idx[12] + 1);
		fprintf(fp, "f %d %d %d %d\n", idx[9] + 1, idx[10] + 1, idx[14] + 1, idx[13] + 1);
		fprintf(fp, "f %d %d %d %d\n", idx[10] + 1, idx[11] + 1, idx[15] + 1, idx[14] + 1);
	}

	fprintf(fp, "o [GregoryPatch]\n");
	for (int i = 0; i < mQuadGregoryPatchIndices.size() / sQuadGregoryPatchSize; i++)
	{
		fprintf(fp, "g GregoryPatchID%04d]\n", i);
		const SizeType* idx = &mQuadGregoryPatchIndices[i * sQuadGregoryPatchSize];
		/*fprintf(fp,
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
				idx[0] + 1, idx[1] + 1, idx[3] + 1, idx[4] + 1, idx[2] + 1,
				idx[5] + 1, idx[6] + 1, idx[8] + 1, idx[9] + 1, idx[7] + 1,
				idx[10] + 1, idx[11] + 1, idx[13] + 1, idx[14] + 1, idx[12] + 1,
				idx[15] + 1, idx[16] + 1, idx[18] + 1, idx[19] + 1, idx[17] + 1,

				idx[1] + 1, idx[7] + 1, idx[9] + 1, idx[3] + 1,
				idx[6] + 1, idx[12] + 1, idx[14] + 1, idx[8] + 1,
				idx[11] + 1, idx[17] + 1, idx[19] + 1, idx[13] + 1,
				idx[16] + 1, idx[2] + 1, idx[4] + 1, idx[18] + 1,

				idx[3] + 1, idx[9] + 1, idx[8] + 1, idx[14] + 1,
				idx[13] + 1, idx[19] + 1, idx[18] + 1, idx[4] + 1);*/
		fprintf(fp,
				// corner
				"f %d %d %d %d %d %d %d %d %d %d %d %d\n"
				"f %d %d %d %d %d %d %d %d\n",
				idx[0] + 1, idx[1] + 1, idx[7] + 1, idx[5] + 1, idx[6] + 1,
				idx[12] + 1, idx[10] + 1, idx[11] + 1, idx[17] + 1, idx[15] + 1,
				idx[16] + 1, idx[2] + 1,
				//
				idx[3] + 1, idx[9] + 1, idx[8] + 1, idx[14] + 1, idx[13] + 1,
				idx[19] + 1, idx[18] + 1, idx[4] + 1);
	}

	fclose(fp);
}

Point3f SubdMesh::faceCenter(SizeType fid) const
{
	auto he = mHDSMesh->heFromFace(fid);
	auto curHE = he;
	Point3f cp;
	SizeType pCount = 0;
	do
	{
		cp += mVerts[curHE->vid];
		pCount++;
		curHE = curHE->next();
	} while (curHE != he);
	if (pCount > 0)
	{
		cp /= pCount;
	}
	return cp;
}

Point3f SubdMesh::edgeCenter(SizeType heid) const
{
	auto he = &mHDSMesh->halfedges[heid];
	return (mVerts[he->vid] + mVerts[he->flip()->vid]) * 0.5f;
}

void SubdMesh::process()
{
	GeomContext context;
	cacheGeometry(context);
	specifyBoundaries(context);
	specifyPatchType(context);
	// Cache vertex point
	computeCornerPoints(context);
	computeEdgePoints(context);

	// Compute result for each patch
	computePatches(context);
}

void SubdMesh::cacheGeometry(GeomContext &context)
{
	// Assume context is cleared before passed in.
	context.mVertexValence.resize(mVerts.size(), 0);
	size_t validFaceCount = mFaceSideCount.size();

	// Cache face centers
	auto &faceCenters = context.mFaceCenter;
	faceCenters.resize(validFaceCount);
	for (size_t i = 0; i < validFaceCount; i++)
	{
		const SizeType faceSide = mFaceSideCount[i];
		const SizeType* vids = &mFaceIndices[mFaceIdOffset[i]];

		for (SizeType j = 0; j < faceSide; j++)
		{
			faceCenters[i] += mVerts[vids[j]];
			context.mVertexValence[vids[j]]++;
		}
		faceCenters[i] /= faceSide;
	}

	// Cache edge middle points
	auto &edgeMiddles = context.mEdgeMiddle;
	edgeMiddles.resize(mHDSMesh->halfedges.size());
	vector<bool> visitedEdges(mHDSMesh->halfedges.size(), false);
	for (auto &he : mHDSMesh->halfedges)
	{
		if (visitedEdges[he.index])
		{
			continue;
		}
		auto hef = he.flip();
		edgeMiddles[he.index] = edgeMiddles[hef->index]
			= (mVerts[he.vid] + mVerts[hef->vid]) * 0.5f;
		visitedEdges[he.index] = visitedEdges[hef->index] = true;
	}
}

void SubdMesh::specifyBoundaries(GeomContext &context)
{

	// Compute vertex valence
	for (size_t i = mFaceSideCount.size(); i < mHDSMesh->faces.size(); i++)
	{
		assert(mHDSMesh->faces[i].isNullFace);
		const HDS::HalfEdge* he = mHDSMesh->heFromFace(i);
		auto curHE = he;
		do
		{
			// If a vertex is on two NullFaces,
			// it's non-manifold and is treated as corner vertex
			if (!context.mBoundaryVertex.insert(curHE->vid).second)
			{
				context.mCornerVertex.insert(curHE->vid);
			}
			curHE = curHE->next();
		} while (curHE != he);
	}
}

void SubdMesh::specifyPatchType(GeomContext &context)
{
	// Assume context is cleared before passed in.
	size_t validFaceCount = mFaceSideCount.size();
	// Process Boundary Vertex
	for (size_t i = validFaceCount; i < mHDSMesh->faces.size(); i++)
	{
		const HDS::HalfEdge* he = mHDSMesh->heFromFace(i);
		auto &curHE = he;
		do
		{
			// If a vertex is on two NullFaces,
			// it's non-manifold and is treated as corner vertex
			if (!context.mBoundaryVertex.insert(curHE->vid).second ||
				context.mVertexValence[curHE->vid] == 1)
			{
				context.mCornerVertex.insert(curHE->vid);
			}
			curHE = curHE->next();
		} while (curHE != he);
	}

	// Find all faces using Gregory Patch
	context.mBadFaceHash.resize(validFaceCount, false);
	for (int i = 0; i < context.mVertexValence.size(); i++)
	{
		if (context.mVertexValence[i] != 4 ||
			context.mBoundaryVertex.find(i) != context.mBoundaryVertex.end())
		{
			const HDS::HalfEdge* he = mHDSMesh->heFromVert(i);
			const HDS::HalfEdge* curHE = he;
			do
			{
				if (curHE->fid < validFaceCount)
				{
					context.mBadFaceHash[curHE->fid] = true;
				}
				curHE = curHE->flip()->next();
			} while (curHE != he);
		}
	}
}

void SubdMesh::computeCornerPoints(const GeomContext &context)
{
	auto &faceCenters = context.mFaceCenter;
	auto &edgeMiddles = context.mEdgeMiddle;
	auto &vertValence = context.mVertexValence;
	auto &cornerVerts = context.mCornerVertex;
	auto &boundaryVerts = context.mBoundaryVertex;

	mPatchVertexBuffer.resize(mVerts.size());
	for (size_t i = 0; i < mVerts.size(); i++)
	{
		const HDS::HalfEdge* he = mHDSMesh->heFromVert(i);
		auto curHE = he;

		const SizeType valence = vertValence[i];
		// Corner vertex
		if (cornerVerts.find(i) != cornerVerts.end())
		{
			mPatchVertexBuffer[i] = mVerts[i];
		}
		// Boundary Vertex
		if (valence == 2 || boundaryVerts.find(i) != boundaryVerts.end())
		{
			// TODO: Valence-2 vertex should be handled properly
			mPatchVertexBuffer[i] += mVerts[i] * 4;
			do
			{
				if (curHE->isBoundary)
				{
					mPatchVertexBuffer[i] += mVerts[curHE->flip()->vid];
				}
				curHE = curHE->rotCCW();
			} while (curHE != he);
			mPatchVertexBuffer[i] /= 6.0f;
		}
		else
		{
			do
			{
				mPatchVertexBuffer[i] +=
					faceCenters[curHE->fid] + edgeMiddles[curHE->index];
				curHE = curHE->rotCCW();
			} while (curHE != he);

			mPatchVertexBuffer[i] =
				mVerts[i] * Gregory::cornerCoef[valence][0] +
				mPatchVertexBuffer[i] * Gregory::cornerCoef[valence][1];
		}
	}
}

void SubdMesh::computeEdgePoints(const GeomContext &context)
{
	// move 2 from q = 2/n * sum(...) to 2/3*lambda q
	const Float cFourOverThree = 4.0 / 3.0;

	size_t initSize = mPatchVertexBuffer.size();
	size_t edgeCount = mHDSMesh->halfedges.size();
	mPatchVertexBuffer.resize(initSize + edgeCount);
	Point3f* edgePoints = &mPatchVertexBuffer[initSize];
	for (auto &he : mHDSMesh->halfedges)
	{
		Point3f q;
		SizeType valence = context.mVertexValence[he.vid];
		const HDS::HalfEdge* curHE = &he;
		for (SizeType i = 0; i < valence; i++)
		{
			q += Gregory::edgeP_coefMi(i, valence) * context.mEdgeMiddle[curHE->index] +
				Gregory::edgeP_coefCi(i, valence) * context.mFaceCenter[curHE->fid];
			curHE = curHE->rotCCW();
		}
		q *= (Gregory::lamda_valence(valence) * cFourOverThree / Float(valence));
		edgePoints[he.index] = q + mPatchVertexBuffer[he.vid];
	}
}

void SubdMesh::computePatches(const GeomContext &context)
{

	for (SizeType fid = 0; fid < mFaceSideCount.size(); fid++)
	{
		const SizeType faceSide = mFaceSideCount[fid];
		//const SizeType* vids = &mFaceIndices[mFaceIdOffset[i]];

		switch (faceSide)
		{
		case 3:
		{
			// generate Triangular Gregory patch
			computeTriGregoryPatch(fid, context);
			break;
		}
		case 4:
		{
			if (context.mBadFaceHash[fid])
			{
				// Generate Quad Gregory Patch
				computeQuadGregoryPatch(fid, context);
			}
			else
			{
				// Generate Bezier Patch
				computeBezierPatch(fid);
			}
			break;
		}
		default:
		{
			break;
		}
		}
	}
}

void SubdMesh::computeBezierPatch(SizeType fid)
{
	std::array<const HDS::HalfEdge*, 4> hes{ mHDSMesh->heFromFace(fid) };
	hes[1] = hes[0]->next();
	hes[2] = hes[1]->next();
	hes[3] = hes[2]->next();

	const SizeType* vids = &mFaceIndices[mFaceIdOffset[fid]];

	assert(hes[0]->vid == vids[0]);

	auto interialPatch = [&](SizeType id1, SizeType id2, SizeType id3, SizeType id4)
	{
		static const Float invNine = 1.0f / 9.0f;
		return (mVerts[id1] * 4 + (mVerts[id2] + mVerts[id4]) * 2 + mVerts[id3]) * invNine;
	};

	SizeType firstFaceId = mPatchVertexBuffer.size();
	mPatchVertexBuffer.resize(firstFaceId + 4);
	mPatchVertexBuffer[firstFaceId] = interialPatch(vids[0], vids[1], vids[2], vids[3]);
	mPatchVertexBuffer[firstFaceId + 1] = interialPatch(vids[1], vids[2], vids[3], vids[0]);
	mPatchVertexBuffer[firstFaceId + 2] = interialPatch(vids[2], vids[3], vids[0], vids[1]);
	mPatchVertexBuffer[firstFaceId + 3] = interialPatch(vids[3], vids[0], vids[1], vids[2]);

	size_t firstPatchId = mBezierPatchIndices.size();
	mBezierPatchIndices.resize(firstPatchId + sQuadBezierPatchSize);

	SizeType* curPatchIdx = &mBezierPatchIndices[firstPatchId];
	curPatchIdx[0] = vids[0];
	curPatchIdx[3] = vids[1];
	curPatchIdx[15] = vids[2];
	curPatchIdx[12] = vids[3];

	SizeType vertexOffset = mVerts.size();
	curPatchIdx[1] = hes[0]->index + vertexOffset;
	curPatchIdx[2] = hes[0]->flip()->index + vertexOffset;
	curPatchIdx[7] = hes[1]->index + vertexOffset;
	curPatchIdx[11] = hes[1]->flip()->index + vertexOffset;
	curPatchIdx[14] = hes[2]->index + vertexOffset;
	curPatchIdx[13] = hes[2]->flip()->index + vertexOffset;
	curPatchIdx[8] = hes[3]->index + vertexOffset;
	curPatchIdx[4] = hes[3]->flip()->index + vertexOffset;

	curPatchIdx[5] = firstFaceId;
	curPatchIdx[6] = firstFaceId + 1;
	curPatchIdx[10] = firstFaceId + 2;
	curPatchIdx[9] = firstFaceId + 3;
}

void SubdMesh::computeQuadGregoryPatch(SizeType fid, const GeomContext &context)
{
	std::array<const HDS::HalfEdge*, 4> hes{ mHDSMesh->heFromFace(fid) };
	hes[1] = hes[0]->next();
	hes[2] = hes[1]->next();
	hes[3] = hes[2]->next();

	const SizeType* vids = &mFaceIndices[mFaceIdOffset[fid]];

	assert(hes[0]->vid == vids[0]);

	SizeType firstFaceId = mPatchVertexBuffer.size();
	mPatchVertexBuffer.resize(firstFaceId + 8);
	Point3f* InteriorFacePoint = &mPatchVertexBuffer[firstFaceId];
	
	const uint32_t firstPatchId = mQuadGregoryPatchIndices.size();
	mQuadGregoryPatchIndices.resize(firstPatchId + sQuadGregoryPatchSize);
	
	SizeType* curPatchIdx = &mQuadGregoryPatchIndices[firstPatchId];
	curPatchIdx[0] = vids[0];
	curPatchIdx[5] = vids[1];
	curPatchIdx[10] = vids[2];
	curPatchIdx[15] = vids[3];

	SizeType vertexOffset = mVerts.size();
	curPatchIdx[1] = hes[0]->index + vertexOffset;
	curPatchIdx[7] = hes[0]->flip()->index + vertexOffset;
	curPatchIdx[6] = hes[1]->index + vertexOffset;
	curPatchIdx[12] = hes[1]->flip()->index + vertexOffset;
	curPatchIdx[11] = hes[2]->index + vertexOffset;
	curPatchIdx[17] = hes[2]->flip()->index + vertexOffset;
	curPatchIdx[16] = hes[3]->index + vertexOffset;
	curPatchIdx[2] = hes[3]->flip()->index + vertexOffset;

	curPatchIdx[3] = firstFaceId;
	curPatchIdx[4] = firstFaceId + 1;
	curPatchIdx[8] = firstFaceId + 2;
	curPatchIdx[9] = firstFaceId + 3;
	curPatchIdx[13] = firstFaceId + 4;
	curPatchIdx[14] = firstFaceId + 5;
	curPatchIdx[18] = firstFaceId + 6;
	curPatchIdx[19] = firstFaceId + 7;

	// Compute Face Points f0+ f0-
	for (SizeType i = 0; i < 4; i++)
	{
		SizeType prevIdx = i == 0 ? 3 : i - 1;
		SizeType nextIdx = i == 3 ? 0 : i + 1;
		// calculate c0,c1,c2
		SizeType curVid = vids[i];
		SizeType prevVid = vids[prevIdx];
		SizeType nextVid = vids[nextIdx];
		auto c0 = Gregory::cosPi(2, context.mVertexValence[curVid]);
		auto c1 = Gregory::cosPi(2, context.mVertexValence[nextVid]);
		auto c2 = Gregory::cosPi(2, context.mVertexValence[prevVid]);

		/*const int dFaceType = 3;// Varies
		const Float inv_d = 1.0 / dFaceType;*/
		const SizeType subPatchSize = 5;
		SizeType subPatchOffset = i * subPatchSize;
		SizeType prevSubPatchOfs = prevIdx * subPatchSize;
		SizeType nextSubPatchOfs = nextIdx * subPatchSize;

		const Point3f &curP0 = mPatchVertexBuffer[curVid];
		const Point3f &e0plus = mPatchVertexBuffer[curPatchIdx[subPatchOffset + 1]];
		const Point3f &e0minus = mPatchVertexBuffer[curPatchIdx[subPatchOffset + 2]];
		const Point3f &e1minus = mPatchVertexBuffer[curPatchIdx[nextSubPatchOfs + 2]];
		const Point3f &e3plus = mPatchVertexBuffer[curPatchIdx[prevSubPatchOfs + 1]];

		const Point3f &m_i = context.mEdgeMiddle[hes[i]->index];
		const Point3f &m_im1 = context.mEdgeMiddle[hes[i]->rotCW()->index];
		const Point3f &m_ip1 = context.mEdgeMiddle[hes[prevIdx]->index];
		const Point3f &m_ip2 = context.mEdgeMiddle[hes[prevIdx]->flip()->prev()->index];

		const Point3f &fc_i = context.mFaceCenter[fid];
		const Point3f &fc_ip1 = context.mFaceCenter[hes[i]->rotCCW()->fid];
		const Point3f &fc_im1 = context.mFaceCenter[hes[i]->flip()->fid];
		// f+
		// get r+
		const Float invThree = 1.0 / 3.0;
		Vector3f r_plus = (m_ip1 - m_im1) * invThree + (fc_i - fc_im1) * 2.0 * invThree;

		InteriorFacePoint[i * 2] = (c1 * curP0
									+ (Gregory::cQuadFaceD - 2 * c0 - c1) * e0plus
									+ 2 * c0 * e1minus
									+ r_plus) * Gregory::cQuadFaceInvD;
		// f-
		// get r-
		Vector3f r_minus = (m_i - m_ip2) * invThree + (fc_i - fc_ip1) * 2.0 * invThree;

		InteriorFacePoint[i * 2 + 1] = (c2 * curP0
										+ (Gregory::cQuadFaceD - 2 * c0 - c2) * e0minus
										+ 2 * c0 * e3plus
										+ r_minus) * Gregory::cQuadFaceInvD;
	}
}

void SubdMesh::computeTriGregoryPatch(SizeType fid, const GeomContext &context)
{
	std::array<const HDS::HalfEdge*, 3> hes{ mHDSMesh->heFromFace(fid) };
	hes[1] = hes[0]->next();
	hes[2] = hes[1]->next();

	const SizeType* vids = &mFaceIndices[mFaceIdOffset[fid]];

	assert(hes[0]->vid == vids[0]);

	SizeType firstFaceId = mPatchVertexBuffer.size();
	mPatchVertexBuffer.resize(firstFaceId + 6);
	Point3f* InteriorFacePoint = &mPatchVertexBuffer[firstFaceId];

	size_t firstPatchId = mTriGregoryPatchIndices.size();
	mTriGregoryPatchIndices.resize(firstPatchId + sQuadBezierPatchSize);

	SizeType* curPatchIdx = &mQuadGregoryPatchIndices[firstPatchId];
	curPatchIdx[0] = vids[0];
	curPatchIdx[5] = vids[1];
	curPatchIdx[10] = vids[2];

	SizeType vertexOffset = mVerts.size();
	curPatchIdx[1] = hes[0]->index + vertexOffset;
	curPatchIdx[7] = hes[0]->flip()->index + vertexOffset;
	curPatchIdx[6] = hes[1]->index + vertexOffset;
	curPatchIdx[12] = hes[1]->flip()->index + vertexOffset;
	curPatchIdx[11] = hes[2]->index + vertexOffset;
	curPatchIdx[2] = hes[2]->flip()->index + vertexOffset;

	curPatchIdx[3] = firstFaceId;
	curPatchIdx[4] = firstFaceId + 1;
	curPatchIdx[8] = firstFaceId + 2;
	curPatchIdx[9] = firstFaceId + 3;
	curPatchIdx[13] = firstFaceId + 4;
	curPatchIdx[14] = firstFaceId + 5;

	// Compute Face Points f0+ f0-
	for (SizeType i = 0; i < 3; i++)
	{
		SizeType prevIdx = i == 0 ? 2 : i - 1;
		SizeType nextIdx = i == 2 ? 0 : i + 1;
		// calculate c0,c1,c2
		SizeType curVid = vids[i];
		SizeType prevVid = vids[prevIdx];
		SizeType nextVid = vids[nextIdx];
		auto c0 = Gregory::cosPi(2, context.mVertexValence[curVid]);
		auto c1 = Gregory::cosPi(2, context.mVertexValence[nextVid]);
		auto c2 = Gregory::cosPi(2, context.mVertexValence[prevVid]);

		const SizeType subPatchSize = 5;
		SizeType subPatchOffset = i * subPatchSize;
		SizeType prevSubPatchOfs = prevIdx * subPatchSize;
		SizeType nextSubPatchOfs = nextIdx * subPatchSize;

		const Point3f &curP0 = mPatchVertexBuffer[curVid];
		const Point3f &e0plus = mPatchVertexBuffer[curPatchIdx[subPatchOffset + 1]];
		const Point3f &e0minus = mPatchVertexBuffer[curPatchIdx[subPatchOffset + 2]];
		const Point3f &e1minus = mPatchVertexBuffer[curPatchIdx[nextSubPatchOfs + 2]];
		const Point3f &e3plus = mPatchVertexBuffer[curPatchIdx[prevSubPatchOfs + 1]];

		const Point3f &m_i = context.mEdgeMiddle[hes[i]->index];
		const Point3f &m_im1 = context.mEdgeMiddle[hes[i]->rotCW()->index];
		const Point3f &m_ip1 = context.mEdgeMiddle[hes[prevIdx]->index];
		const Point3f &m_ip2 = context.mEdgeMiddle[hes[prevIdx]->flip()->prev()->index];

		const Point3f &fc_i = context.mFaceCenter[fid];
		const Point3f &fc_ip1 = context.mFaceCenter[hes[i]->rotCCW()->fid];
		const Point3f &fc_im1 = context.mFaceCenter[hes[i]->flip()->fid];
		// f+
		// get r+
		const Float invThree = 1.0 / 3.0;
		Vector3f r_plus = (m_ip1 - m_im1) * invThree + (fc_i - fc_im1) * 2.0 * invThree;

		InteriorFacePoint[i * 2] = (c1 * curP0
									+ (Gregory::cTriFaceD - 2 * c0 - c1) * e0plus
									+ 2 * c0 * e1minus
									+ r_plus) * Gregory::cTriFaceInvD;
		// f-
		// get r-
		Vector3f r_minus = (m_i - m_ip2) * invThree + (fc_i - fc_ip1) * 2.0 * invThree;

		InteriorFacePoint[i * 2 + 1] = (c2 * curP0
										+ (Gregory::cTriFaceD - 2 * c0 - c2) * e0minus
										+ 2 * c0 * e3plus
										+ r_minus) * Gregory::cTriFaceInvD;
	}
}

#if 0

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
						   float s, float t) {
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
					mGregoryPatch[0],
					mGregoryPatch[5],
					mGregoryPatch[10],
					mGregoryPatch[15],
					u, v
				));
			}
			else
			{
				// Sample on grid
				// Face Points
				Point3f f0 = (u * mGregoryPatch[3]
							  + v * mGregoryPatch[4])
					/ (u + v);
				Point3f f1 = ((1 - u) * mGregoryPatch[9]
							  + v  * mGregoryPatch[8])
					/ (1 - u + v);
				Point3f f2 = ((1 - u) * mGregoryPatch[13]
							  + (1 - v) * mGregoryPatch[14])
					/ (2 - u - v);
				Point3f f3 = (u * mGregoryPatch[19]
							  + (1 - v) * mGregoryPatch[18])
					/ (1 + u - v);

				Point3f p[9];
				// 6 -- 7 -- 8
				// |    |    |
				// 3 -- 4 -- 5
				// |    |    |
				// 0 -- 1 -- 2
				p[0] = bilinear_pos(mGregoryPatch[0],
									mGregoryPatch[1],
									f0,
									mGregoryPatch[2],
									u, v);
				p[1] = bilinear_pos(mGregoryPatch[1],
									mGregoryPatch[7],
									f1,
									f0,
									u, v);
				p[2] = bilinear_pos(mGregoryPatch[7],
									mGregoryPatch[5],
									mGregoryPatch[6],
									f1,
									u, v);

				p[3] = bilinear_pos(mGregoryPatch[2],
									f0,
									f3,
									mGregoryPatch[16],
									u, v);
				p[4] = bilinear_pos(f0, f1, f2, f3,
									u, v);
				p[5] = bilinear_pos(f1,
									mGregoryPatch[6],
									mGregoryPatch[12],
									f2,
									u, v);

				p[6] = bilinear_pos(mGregoryPatch[16],
									f3,
									mGregoryPatch[17],
									mGregoryPatch[15],
									u, v);
				p[7] = bilinear_pos(f3,
									f2,
									mGregoryPatch[11],
									mGregoryPatch[17],
									u, v);
				p[8] = bilinear_pos(f2,
									mGregoryPatch[12],
									mGregoryPatch[10],
									mGregoryPatch[11],
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
#endif

}
