#version 450
layout(location = 0) in vec3 vp;

uniform mat4 view_matrix;

out vec3 position;

void main()
{
	gl_Position = view_matrix * vec4(vp, 1.0);
	position = gl_Position.xyz;
}