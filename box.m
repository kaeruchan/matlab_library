function [x,y] = box(r1,r2,s)
% box muller
  x = s * sqrt(-2.0*log(r1)) * cos(2.0*PI*r2);
  y = s * sqrt(-2.0*log(r1)) * sin(2.0*PI*r2);
 end
