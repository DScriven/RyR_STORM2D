#version 430 core
layout (location=0) in vec2 vpos;
layout (location=1) in float logDensity;
out vec4 pcolor;
uniform float PointSize;
uniform mat4 transform;
uniform vec4 DensityColor[8];
void main(void)
{
    int cno = 0;
    if(logDensity < -3.0)
        cno = 1;
    else if (logDensity >= -3.0 && logDensity < -2.5)
        cno = 1;
    else if (logDensity >= -2.5 && logDensity < -2)
        cno = 2;
    else if (logDensity >= -2 && logDensity < -1.75)
        cno = 3;
    else if (logDensity >= -1.75 && logDensity < -1.5)
        cno = 4;
    else if (logDensity >= -1.5 && logDensity < -1.25)
        cno = 5;
    else if (logDensity >= -1.25 && logDensity < -1.12)
        cno = 6;
    else if (logDensity >= -1.12)
        cno = 7;

    pcolor = vec4(DensityColor[cno]);
    gl_PointSize = PointSize;
    gl_Position = transform*vec4(vpos,0.0,1.0);
}

