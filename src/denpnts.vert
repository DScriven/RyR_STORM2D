#version 430 core
layout (location=0) in vec2 vpos;
layout (location=1) in float logDensity;
out vec4 pcolor;
uniform int AboveThresholdColourIndex;
uniform int BelowThresholdColourIndex;
uniform float BlinkSize;
uniform float denThreshold;
uniform mat4 transform;
uniform vec4 DensityColor[8];
void main(void)
{
    int cno = AboveThresholdColourIndex;

    if(logDensity < denThreshold)
        cno = BelowThresholdColourIndex;

    pcolor = vec4(DensityColor[cno]);
    gl_PointSize = BlinkSize;
    gl_Position = transform*vec4(vpos,0.0,1.0);
}

