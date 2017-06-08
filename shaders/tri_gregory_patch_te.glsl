#version 450
layout(triangles, equal_spacing) in;

uniform mat4 proj_matrix;

patch in vec3 patchPts[15];
patch in vec3 edge[12];

out vec3 pos_eye, norm_eye;

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


void main()
{
	float u = gl_TessCoord.x;// (n_layer - 1) * t;
	float v = gl_TessCoord.y;
	float w = gl_TessCoord.z;

	// TODO need to be re-written
	
	if ((u == 0 || u == 1) && (v == 0 || v == 1) && (w == 0 || w == 1))
	{
		int idx = u == 1 ? 0 : (v == 1 ? 5 : 10);
		/*if (u == 1 && v == 0 && w == 0)
			idx = 0;
		if (u == 0 && v == 1 && w == 0)
			idx = 5;
		if (u == 0 && v == 0 && w == 1)
			idx = 10;*/
		pos_eye = patchPts[idx];
		vec3 dpdu = patchPts[idx+1] - patchPts[idx];
		vec3 dpdv = patchPts[idx+2] - patchPts[idx];
		norm_eye = normalize(cross(dpdu, dpdv));
		gl_Position = proj_matrix * vec4(pos_eye, 1.0f);
	}
	else
	{
		// Sample on grid
		// Face Points
		vec3 f0 = (v * patchPts[ 3] + w * patchPts[ 4]) / (v + w);
		vec3 f1 = (w * patchPts[ 9] + u * patchPts[ 8]) / (w + u);
		vec3 f2 = (u * patchPts[13] + v * patchPts[14])	/ (u + v);

		vec3 p[10];
		//          09
		//         /  \
		//       05 -- 08
		//      /  \  /  \
		//    02 -- 04 -- 07
		//   /  \  /  \  /  \
		// 00 -- 01 -- 03 -- 06
		p[0] = u * edge[ 0].xyz + v * edge[ 1].xyz + w * edge[11].xyz;
		p[1] = u * edge[ 1].xyz + v * edge[ 2].xyz + w * f0;
		p[2] = u * edge[11].xyz + v * f0           + w * edge[10].xyz;
		p[3] = u * edge[ 2].xyz + v * edge[ 3].xyz + w * f1;
		p[4] = u * f0           + v * f1           + w * f2;
		p[5] = u * edge[10].xyz + v * f2           + w * edge[ 9].xyz;
		p[6] = u * edge[ 3].xyz + v * edge[ 4].xyz + w * edge[ 5].xyz;
		p[7] = u * f1           + v * edge[ 5].xyz + w * edge[ 6].xyz;
		p[8] = u * f2           + v * edge[ 6].xyz + w * edge[ 7].xyz;
		p[9] = u * edge[ 9].xyz + v * edge[ 7].xyz + w * edge[ 8].xyz;

		int idx[4] = { 0, 1, 3, 6 };
		for (int i = 3; i > 1; i--)
		{
			for (int j = 0; j < i; j++)
			{
				for (int k = 0; k <= j; k++)
				{
					p[idx[j] + k] = u * p[idx[j] + k] + v * p[idx[j+1]+k] + w * p[idx[j+1]+k + 1];
				}
			}
		}
		
		pos_eye = u * p[0] + v * p[1] + w * p[2];
		vec3 dpdu = p[1] - p[0];
		vec3 dpdv = p[2] - p[0];

		norm_eye = normalize(cross(dpdu, dpdv));

		gl_Position = proj_matrix * vec4(pos_eye, 1.0f);
	}
}