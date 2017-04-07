#define MAXFLOAT ((float)3.40282346638528860e+38)

int
finite(double x)
{
  if (x < MAXFLOAT && x > - MAXFLOAT)  return(1);
  else return(0);
}


double 
vrint(double x)
{
  double y, iy;
  int i;

  i  = (int) x;
  iy = (double) i;

  if (x >= 0) {
     if (iy + 0.5 > x) y = iy;
     else y = iy + 1.0;
  }
  else {
    if (iy - 0.5 < x) y = iy;
    else y = iy - 1.0;
  }
  return(y);
}

 
float 
vrintf(float x)
{
  float y, iy;
  int i;

  i  = (int) x;
  iy = (float) i;

  if (x >= 0) {
     if (iy + 0.5 > x) y = iy;
     else y = iy + 1.0;
  }
  else {
    if (iy - 0.5 < x) y = iy;
    else y = iy - 1.0;
  }
  return(y);
}
