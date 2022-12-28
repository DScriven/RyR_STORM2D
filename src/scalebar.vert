in vec2 vpos;
out vec4 pcolor;
uniform mat4 transform;
uniform vec4 ScaleBarColour;

void main(void)
{
    pcolor = ScaleBarColour;
    gl_Position = transform * vec4(vpos,0.0,1.0);
}
