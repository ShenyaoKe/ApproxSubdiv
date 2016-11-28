#version 430

layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;

in vec3 pos_eye_te[];
out vec3 pos_eye, norm_eye;

void main()
{
	// Transform each vertex into viewport space
	vec4 p0 = gl_in[0].gl_Position;
	vec4 p1 = gl_in[1].gl_Position;
	vec4 p2 = gl_in[2].gl_Position;

	vec4 v1 = normalize(gl_in[1].gl_Position - gl_in[0].gl_Position);
	vec4 v2 = normalize(gl_in[2].gl_Position - gl_in[0].gl_Position);

	//vec3 norm = normalize(cross(v1.xyz, v2.xyz));
	//norm_eye = normalize(view_matrix * vec4(norm, 0)).xyz;
	norm_eye = normalize(cross(v1.xyz, v2.xyz));
	//float d[4];
	// Send the triangle along with the edge distances
	gl_PrimitiveID = gl_PrimitiveIDIn; 
	
	gl_Position = p0;
	pos_eye = pos_eye_te[0];
	EmitVertex();
	gl_Position = p1;
	pos_eye = pos_eye_te[1];
	EmitVertex();
	gl_Position = p2;
	pos_eye = pos_eye_te[2];
	EmitVertex();
}