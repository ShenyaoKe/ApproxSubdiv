#version 450
layout(quads, equal_spacing) in;

uniform mat4 proj_matrix;

out vec3 pos_eye, norm_eye;

vec3 bilinear_pos(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float s, float t)
{
	// id3 -- id2
	//  |      |
	// id0 -- id1
	return mix(mix(p0, p1, s), mix(p3, p2, s), t);
}

// Gregory Patch Layout
// 15-----17-----11-----10
// |      |      |      |
// |      19     13     |
// 16--18          14--12
// |       F3  F2       |
// |       F0  F1       |
// 02--04          08--06
// |      03     09     |
// |      |      |      |
// 00-----01-----07-----05
void main()
{
	float u = gl_TessCoord.x;// (n_layer - 1) * t;
	float v = gl_TessCoord.y;
	/*if ((u == 0 || u == 1) && (v == 0 || v == 1))
	{
		pos_eye = bilinear_pos(
									gl_in[ 0].gl_Position.xyz,
									gl_in[ 5].gl_Position.xyz,
									gl_in[10].gl_Position.xyz,
									gl_in[15].gl_Position.xyz,
									u, v
									);
		gl_Position = proj_matrix * vec4(pos_eye, 1.0f);
	}*/
	if (u == 0 && v == 0)
	{
		pos_eye = gl_in[ 0].gl_Position.xyz;
		vec4 dpdu = gl_in[ 3].gl_Position - gl_in[ 0].gl_Position;
		vec4 dpdv = gl_in[ 2].gl_Position - gl_in[ 0].gl_Position;
		norm_eye = normalize(cross(dpdu.xyz, dpdv.xyz));
		gl_Position = proj_matrix * gl_in[0].gl_Position;
	}
	else if (u == 1 && v == 0)
	{
		pos_eye = gl_in[ 5].gl_Position.xyz;
		vec4 dpdu = gl_in[ 6].gl_Position - gl_in[ 5].gl_Position;
		vec4 dpdv = gl_in[ 7].gl_Position - gl_in[ 5].gl_Position;
		norm_eye = normalize(cross(dpdu.xyz, dpdv.xyz));
		gl_Position = proj_matrix * gl_in[5].gl_Position;
	}
	else if (u == 1 && v == 1)
	{
		pos_eye = gl_in[10].gl_Position.xyz;
		vec4 dpdu = gl_in[11].gl_Position - gl_in[10].gl_Position;
		vec4 dpdv = gl_in[12].gl_Position - gl_in[10].gl_Position;
		norm_eye = normalize(cross(dpdu.xyz, dpdv.xyz));
		gl_Position = proj_matrix * gl_in[10].gl_Position;
	}
	else if (u == 0 && v == 1)
	{
		pos_eye = gl_in[15].gl_Position.xyz;
		vec4 dpdu = gl_in[16].gl_Position - gl_in[15].gl_Position;
		vec4 dpdv = gl_in[17].gl_Position - gl_in[15].gl_Position;
		norm_eye = normalize(cross(dpdu.xyz, dpdv.xyz));
		gl_Position = proj_matrix * gl_in[15].gl_Position;
	}
	else
	{
		// Sample on grid
		// Face Points
		vec3 f0 = (      u * gl_in[ 3].gl_Position
				+        v * gl_in[ 4].gl_Position).xyz
				/ (u + v);
		vec3 f1 = ((1 - u) * gl_in[ 9].gl_Position
				+       v  * gl_in[ 8].gl_Position).xyz
				/ (1 - u + v);
		vec3 f2 = ((1 - u) * gl_in[13].gl_Position
				+  (1 - v) * gl_in[14].gl_Position).xyz
				/ (2 - u - v);
		vec3 f3 = (      u * gl_in[19].gl_Position
				+  (1 - v) * gl_in[18].gl_Position).xyz
				/ (1 + u - v);

		vec3 p[9];
		// 6 -- 7 -- 8
		// |    |    |
		// 3 -- 4 -- 5
		// |    |    |
		// 0 -- 1 -- 2
		p[0] = bilinear_pos(gl_in[0].gl_Position.xyz,
						gl_in[1].gl_Position.xyz,
						f0,
						gl_in[2].gl_Position.xyz,
						u, v);
		p[1] = bilinear_pos(gl_in[1].gl_Position.xyz,
						gl_in[7].gl_Position.xyz,
						f1,
						f0,
						u, v);
		p[2] = bilinear_pos(gl_in[7].gl_Position.xyz,
						gl_in[5].gl_Position.xyz,
						gl_in[6].gl_Position.xyz,
						f1,
						u, v);

		p[3] = bilinear_pos(gl_in[2].gl_Position.xyz,
						f0,
						f3,
						gl_in[16].gl_Position.xyz,
						u, v);
		p[4] = bilinear_pos(f0, f1, f2, f3,
						u, v);
		p[5] = bilinear_pos(f1,
						gl_in[6].gl_Position.xyz,
						gl_in[12].gl_Position.xyz,
						f2,
						u, v);

		p[6] = bilinear_pos(gl_in[16].gl_Position.xyz,
						f3,
						gl_in[17].gl_Position.xyz,
						gl_in[15].gl_Position.xyz,
						u, v);
		p[7] = bilinear_pos(f3,
						f2,
						gl_in[11].gl_Position.xyz,
						gl_in[17].gl_Position.xyz,
						u, v);
		p[8] = bilinear_pos(f2,
						gl_in[12].gl_Position.xyz,
						gl_in[10].gl_Position.xyz,
						gl_in[11].gl_Position.xyz,
						u, v);
	
		// 2 -- 3
		// |    |
		// 0 -- 1
		p[0] = bilinear_pos(p[0], p[1], p[4], p[3], u, v);
		p[1] = bilinear_pos(p[1], p[2], p[5], p[4], u, v);
		p[2] = bilinear_pos(p[3], p[4], p[7], p[6], u, v);
		p[3] = bilinear_pos(p[4], p[5], p[8], p[7], u, v);

		pos_eye = bilinear_pos(p[0], p[1], p[3], p[2], u, v);

		vec3 dpdu = mix(p[1] - p[0], p[3] - p[2], v);
		vec3 dpdv = mix(p[2] - p[0], p[3] - p[1], u);
		norm_eye = normalize(cross(dpdu, dpdv));

		gl_Position = proj_matrix * vec4(pos_eye, 1.0f);
	}
}