#version 430

layout(lines_adjacency) in;
layout(triangle_strip, max_vertices = 4) out;

//out vec3 norm;
//out vec2 TexCoord;
out vec3 pos_eye, norm_eye;
out float cosTheta;

uniform mat4 view_matrix, proj_matrix;
//uniform vec3 lpos_world = vec3(0.0, 0, 0);

void main()
{
	mat4 mvp = proj_matrix * view_matrix;
	// Transform each vertex into viewport space
	vec4 p0 = mvp * gl_in[0].gl_Position;
	vec4 p1 = mvp * gl_in[1].gl_Position;
	vec4 p2 = mvp * gl_in[2].gl_Position;
	vec4 p3 = mvp * gl_in[3].gl_Position;

	vec4 v1 = normalize(gl_in[2].gl_Position - gl_in[0].gl_Position);
	vec4 v2 = normalize(gl_in[3].gl_Position - gl_in[1].gl_Position);

	vec3 norm = normalize(cross(v1.xyz, v2.xyz));
	norm_eye = normalize((view_matrix * vec4(norm, 0)).xyz);

	// Send the triangle along with the edge distances
	gl_PrimitiveID = gl_PrimitiveIDIn; 
	
	gl_Position = p0;
	pos_eye = (view_matrix * gl_in[0].gl_Position).xyz;
	//cosTheta = dot(normalize(lpos_world - pos_eye), norm_eye);
	EmitVertex();

	gl_Position = p1;
	pos_eye = (view_matrix * gl_in[1].gl_Position).xyz;
	//cosTheta = dot(normalize(lpos_world - pos_eye), norm_eye);
	EmitVertex();

	gl_Position = p3;
	pos_eye = (view_matrix * gl_in[3].gl_Position).xyz;
	//cosTheta = dot(normalize(lpos_world - pos_eye), norm_eye);
	EmitVertex();

	gl_Position = p2;
	pos_eye = (view_matrix * gl_in[2].gl_Position).xyz;
	//cosTheta = dot(normalize(lpos_world - pos_eye), norm_eye);
	EmitVertex();
}