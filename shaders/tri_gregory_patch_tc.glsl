#version 450
layout(vertices = 1) out;

uniform float segments = 16.0f;

in vec3 position[];
patch out vec3 patchPts[15];
patch out vec3 edge[12];

// Triangular Gregory Patch Layout
//               10                                    08
//             /    \                                /    \
//           11      12                            09      07
//          /  \    /  \                          /          \
//        /    13  14    \        ==>           10            06
//      02--04        08--06                   /                \
//     /      03    09      \                11                  05
//   /         |    |         \             /                      \
//  00--------01----07--------05          00----01-----02-----03----04

vec3 elevate1(vec3 p1, vec3 p2)
{
	return 0.25f * p1 + 0.75f * p2;
}
vec3 elevate2(vec3 p1, vec3 p2)
{
	return 0.5f * p1 + 0.5f * p2;
}
vec3 elevate3(vec3 p1, vec3 p2)
{
	return 0.75f * p1 + 0.25f * p2;
}

void main()
{
	for (int i = 0; i < 15; i++)
	{
		patchPts[i] = position[i];
	}
	edge[0] = position[0];
	edge[1] = elevate1(position[0], position[1]);
	edge[2] = elevate2(position[1], position[7]);
	edge[3] = elevate3(position[7], position[5]);
	edge[4] = position[5];
	edge[5] = elevate1(position[5], position[6]);
	edge[6] = elevate2(position[6], position[12]);
	edge[7] = elevate3(position[12], position[10]);
	edge[8] = position[10];
	edge[9] = elevate1(position[10], position[11]);
	edge[10] = elevate2(position[11], position[2]);
	edge[11] = elevate3(position[2], position[0]);
	

	// Define the tessellation levels
	gl_TessLevelInner[0] = segments;
	gl_TessLevelOuter[0] = segments;
	gl_TessLevelOuter[1] = segments;
	gl_TessLevelOuter[2] = segments;

}