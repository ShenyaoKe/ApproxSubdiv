#version 430
/*in vec3 normal;*/
uniform unsigned int rootColor, tipColor, scatterColor;// = 0xFFFFFF;
in vec3 uvw;
out vec4 frag_color;
void main ()
{
	vec3 root_color = vec3(
		float((rootColor >> 16) & 0xFF) / 255.0,
		float((rootColor >> 8) & 0xFF) / 255.0,
		float(rootColor & 0xFF) / 255.0);
	vec3 tip_color = vec3(
		float((tipColor >> 16) & 0xFF) / 255.0,
		float((tipColor >> 8) & 0xFF) / 255.0,
		float(tipColor & 0xFF) / 255.0);
	float alpha = (uvw.z - 0.8) / 0.2;
	alpha = alpha > 1 ? 1 : (alpha < 0 ? 0 : alpha);
	frag_color = vec4(mix(root_color, tip_color, (uvw.z * uvw.z)), 1-alpha);
}