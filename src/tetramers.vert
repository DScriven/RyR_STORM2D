#version 430 core
layout (location=0) in vec2 vpos;
layout (location=1) in vec3 linecolour;
out vec4 pcolor;
uniform float transparency;
uniform mat4 transform;
void main(void)
{
    pcolor = vec4(linecolour,transparency);
    gl_Position = transform*vec4(vpos,0.0,1.0);
}

