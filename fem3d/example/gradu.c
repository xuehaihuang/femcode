  t2 = pii*x;
  t3 = pii*y;
  t4 = pii*z;
  t5 = pii*pii;
  t6 = cos(t2);
  t7 = cos(t3);
  t8 = cos(t4);
  t9 = sin(t2);
  t10 = sin(t3);
  t11 = sin(t4);
  t12 = t9*t9;
  t13 = t10*t10;
  t14 = t11*t11;
  t15 = t5*t6*t7*t9*t10*t14*4.0;
  A0[0][0] = t15;
  A0[0][1] = t5*t12*t14*cos(t3*2.0)*2.0;
  A0[0][2] = t5*t7*t8*t10*t11*t12*4.0;
  A0[1][0] = t5*t13*t14*cos(t2*2.0)*-2.0;
  A0[1][1] = -t15;
  A0[1][2] = t5*t6*t8*t9*t11*t13*-4.0;