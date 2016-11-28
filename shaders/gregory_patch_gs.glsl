#version 450

layout(triangles) in;
layout(triangle_strip, max_vertices = 4) out;
uniform mat4 proj_matrix;

out vec3 pos_eye, norm_eye;

void main()
{
	// Transform each vertex into viewport space
	vec4 p0 = gl_in[0].gl_Position;
	vec4 p1 = gl_in[1].gl_Position;
	vec4 p2 = gl_in[2].gl_Position;

	vec4 v1 = p1 - p0;
	vec4 v2 = p2 - p0;

	//vec3 norm = normalize(cross(v1.xyz, v2.xyz));
	//norm_eye = normalize(view_matrix * vec4(norm, 0)).xyz;
	norm_eye = normalize(cross(v1.xyz, v2.xyz));
	//float d[4];
	// Send the triangle along with the edge distances
	gl_PrimitiveID = gl_PrimitiveIDIn; 
	
	gl_Position = proj_matrix * p0;
	pos_eye = p0.xyz;
	EmitVertex();

	gl_Position = proj_matrix * p1;
	pos_eye = p1.xyz;
	EmitVertex();

	gl_Position = proj_matrix * p2;
	pos_eye = p2.xyz;
	EmitVertex();
}