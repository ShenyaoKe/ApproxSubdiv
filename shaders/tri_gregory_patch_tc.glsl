#version 450
layout(vertices = 15) out;

uniform float segments = 16.0f;
uniform int nStrips = 1;// Number of Strips
//patch out int vertCount;

void main()
{
	//vertCount = gl_PatchVerticesIn;
	// Pass along the vertex position unmodified
	gl_out[gl_InvocationID].gl_Position = gl_in[gl_InvocationID].gl_Position;
	// Define the tessellation levels
	gl_TessLevelInner[0] = segments;
	gl_TessLevelInner[1] = segments;
	gl_TessLevelOuter[0] = segments;
	gl_TessLevelOuter[1] = segments;
	gl_TessLevelOuter[2] = segments;
	gl_TessLevelOuter[3] = segments;
}