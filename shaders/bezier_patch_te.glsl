#version 450
layout(quads) in;

uniform mat4 proj_matrix;

out vec3 pos_eye, norm_eye;

vec3 bezier(int id0, int id1, int id2, int id3, float t)
{
	vec3 p[3];// Control points

	p[0] = mix(gl_in[id0].gl_Position, gl_in[id1].gl_Position, t).xyz;
	p[1] = mix(gl_in[id1].gl_Position, gl_in[id2].gl_Position, t).xyz;
	p[2] = mix(gl_in[id2].gl_Position, gl_in[id3].gl_Position, t).xyz;
	
	p[0] = mix(p[0], p[1], t);
	p[1] = mix(p[1], p[2], t);
	
	return mix(p[0], p[1], t);
}
vec3 bilinear_ori(int id0, int id1, int id2, int id3, float s, float t)
{
	// id3 -- id2
	//  |      |
	// id0 -- id1
	return mix( mix(gl_in[id0].gl_Position, gl_in[id1].gl_Position, s),
				mix(gl_in[id3].gl_Position, gl_in[id2].gl_Position, s),
				t).xyz;
}



vec3 bilinear(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float s, float t)
{
	// id3 -- id2
	//  |      |
	// id0 -- id1
	return mix(mix(p0, p1, s), mix(p3, p2, s), t);
}

void main()
{
	float u = gl_TessCoord.x;// (n_layer - 1) * t;
	float v = gl_TessCoord.y;
	// Sample on grid
#if 0
	vec3 p[4];
	p[0] = bezier(0, 1, 2, 3, u);
	p[1] = bezier(4, 5, 6, 7, u);
	p[2] = bezier(8, 9, 10, 11, u);
	p[3] = bezier(12, 13, 14, 15, u);
	
	p[0] = mix(p[0], p[1], v);
	p[1] = mix(p[1], p[2], v);
	p[2] = mix(p[2], p[3], v);
	
	
	p[0] = mix(p[0], p[1], v);
	p[1] = mix(p[1], p[2], v);
	pos_eye = mix(p[0], p[1], v);

	vec3 dpdv = p[1] - p[0];
	vec3 dpdu = (gl_in[1].gl_Position - gl_in[0].gl_Position).xyz;
#else
	vec3 p[9];
	// 6 -- 7 -- 8
	// |    |    |
	// 3 -- 4 -- 5
	// |    |    |
	// 0 -- 1 -- 2
	p[0] = bilinear_ori(0, 1, 5, 4, u, v);
	p[1] = bilinear_ori(1, 2, 6, 5, u, v);
	p[2] = bilinear_ori(2, 3, 7, 6, u, v);

	p[3] = bilinear_ori(4, 5, 9, 8, u, v);
	p[4] = bilinear_ori(5, 6, 10, 9, u, v);
	p[5] = bilinear_ori(6, 7, 11, 10, u, v);

	p[6] = bilinear_ori(8, 9, 13, 12, u, v);
	p[7] = bilinear_ori(9, 10, 14, 13, u, v);
	p[8] = bilinear_ori(10, 11, 15, 14, u, v);
	
	// 2 -- 3
	// |    |
	// 0 -- 1
	p[0] = bilinear(p[0], p[1], p[4], p[3], u, v);
	p[1] = bilinear(p[1], p[2], p[5], p[4], u, v);
	p[2] = bilinear(p[3], p[4], p[7], p[6], u, v);
	p[3] = bilinear(p[4], p[5], p[8], p[7], u, v);

	pos_eye = bilinear(p[0], p[1], p[3], p[2], u, v);
	
	vec3 dpdu = 2 * mix(p[1] - p[0], p[3] - p[2], v);
	vec3 dpdv = 2 * mix(p[2] - p[0], p[3] - p[1], u);
#endif
	

	norm_eye = normalize(cross(dpdu, dpdv));

	gl_Position = proj_matrix * vec4(pos_eye, 1.0f);
}