#version 440
layout(isolines) in;

uniform mat4 view_matrix, proj_matrix;

uniform int nLayer;
patch in int vertCount;
patch in int gridSample;
out vec3 uvw;

vec4 p[4];// Control points

vec4 posAtLayer(int layer, float s, float t)
{
	int idx = layer * 4;

	return mix(
		mix(gl_in[idx].gl_Position, gl_in[idx + 1].gl_Position, s),
		mix(gl_in[idx + 3].gl_Position, gl_in[idx + 2].gl_Position, s),
		t);
}
float rand(vec2 co)
{
	return (fract(sin(dot(co.xy, vec2(12.9898, 78.233))) * 43758.5453) + 1) * 0.5;
}
vec4 blossomCatmullRom(float t)
{
	// Catmull-Rom Spline
	t = t + 1; // 1: degree - 1 = 2 - 1

	// Lagrange
	// Degree: 2
	for (int j = 0; j < 2; j++)
	{
		for (int i = 0; i < 3 - j; i++)
		{
			p[i] = p[i] * (i + j + 1 - t) / float(j + 1)
				+ p[i + 1] * (t - i) / float(j + 1);
		}
	}
	//de Boor
	return p[0] * (2 - t) + p[1] * (t - 1);
}
void main()
{

	float u = gl_TessCoord.x;// (n_layer - 1) * t;
	float v = gl_TessCoord.y;
	// Sample on grid
	int gridSampleSq = gridSample * gridSample;
	int planar_index = int(v * gridSampleSq);
	float si = float(planar_index % gridSample) + rand(vec2(v, v));
	float ti = float(planar_index / gridSample) + rand(vec2(v, v));
	float s = si / float(gridSample);
	float t = ti / float(gridSample);

	uvw = vec3(s, t, u);

	// Catmull-Rom Blossom
	float cvT = u * float(nLayer - 1);
	int iLayer = int(cvT);
	int jLayer = iLayer + 1;
	int prevLayer = iLayer - 1;
	int nextLayer = jLayer + 1;
	p[1] = posAtLayer(iLayer, s, t);
	p[2] = posAtLayer(jLayer, s, t);
	p[0] = prevLayer < 0 ? p[1] : posAtLayer(prevLayer, s, t);
	p[3] = nextLayer >= nLayer ? p[2] : posAtLayer(nextLayer, s, t);
	
	cvT -= iLayer;
	
	gl_Position = proj_matrix * view_matrix * blossomCatmullRom(cvT);
}