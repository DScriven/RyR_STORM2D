in float density
vec3 dcolor
void main(void)
{
      if (density < -3)
        dcolor = (1.0, 1.0, 1.0);
      elseif (density >= -3 && density < -2.5)
        dcolor = (0.0, 0.0, 1.0);
      elseif (density >= -2.5 && density < -2)
        dcolor = (0.0, 1.0, 1.0);
      elseif (density >= -2 && density < -1.75)
        dcolor = (0.0, 1.0, 0.0);
      elseif (density >= -1.75 && density < -1.5)
        dcolor = [1.0 0.65 0.0];
      elseif (density >= -1.5 && density < -1.25)
        dcolor = [1.0 0.0 0.0];
      elseif (density >= -1.25 && density < -1)
        dcolor = [1.0 0.0 1.0];
      elseif (density >= -1)
        dcolor = (0.0, 0.0, 0.0);
      end
      gl_FragColor = dcolor*glPosition;
}
