  t2 = pii*x;
  t3 = pii*y;
  t4 = pii*z;
  t5 = pii*pii;
  t6 = cos(t4);
  t7 = sin(t2);
  t8 = sin(t3);
  t9 = sin(t4);
  t10 = t7*t7;
  t11 = t8*t8;
  A0[0][0] = t5*t6*t7*t9*t11*cos(t2)*4.0;
  A0[1][0] = t5*t6*t8*t9*t10*cos(t3)*4.0;
  A0[2][0] = t5*(t9*t9)*(t10+t11-t10*t11*4.0)*-2.0;
