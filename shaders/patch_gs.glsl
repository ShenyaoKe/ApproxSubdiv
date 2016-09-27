#version 450
// extension should be core -- commented out
// replaced gl_VerticesIn with gl_in.length ()
//#extension GL_EXT_geometry_shader4 : enable

layout(lines) in;
// convert to points, line_strip, or triangle_strip
layout(triangle_strip, max_vertices = 5) out;
//layout(line_strip) out;
uniform mat4 proj_matrix;
in vec2 uv_te[];
//in vec4 T_te[];
//in vec4 B_te[];

out vec2 uv_gs;
uniform float pointsize = 0.02;

void main ()
{
	
	//vec4 cur_pos = gl_in[0].gl_Position;
	vec4 T = normalize(gl_in[1].gl_Position - gl_in[0].gl_Position);
	vec4 N = vec4(0, 0, 1, 0);
	vec4 B = vec4(cross(T.xyz, N.xyz), 0);

	gl_Position = proj_matrix * (gl_in[0].gl_Position - B_te[0] * pointsize);
	uv_gs = uv_te[0];
	EmitVertex();
	gl_Position = proj_matrix * (gl_in[1].gl_Position - B_te[1] * pointsize);
	uv_gs = uv_te[1];
	EmitVertex();
	gl_Position = proj_matrix * (gl_in[0].gl_Position + B_te[0] * pointsize);
	uv_gs = vec2(uv_te[0].x, 1.0);
	EmitVertex();
	gl_Position = proj_matrix * (gl_in[1].gl_Position + B_te[1] * pointsize);
	uv_gs = vec2(uv_te[1].x, 1.0);
	EmitVertex();
}
