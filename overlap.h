/* Author: Stefan Gottschalk (extracted form his PhD thesis)
 */
#ifndef OVERLAP_H
#define OVERLAP_H

#define myfabs(x) ((x < 0) ? -x : x)

template< typename NT >
int
obb_disjoint(const NT B[3][3], const NT T[3], const NT a[3], const NT b[3])
{
  /*register*/ NT t, s;
  /*register*/ int r;
  NT Bf[3][3];
  const NT reps = 1e-6;
  
  // Bf = fabs(B)
  Bf[0][0] = myfabs(B[0][0]);  Bf[0][0] += reps;
  Bf[0][1] = myfabs(B[0][1]);  Bf[0][1] += reps;
  Bf[0][2] = myfabs(B[0][2]);  Bf[0][2] += reps;
  Bf[1][0] = myfabs(B[1][0]);  Bf[1][0] += reps;
  Bf[1][1] = myfabs(B[1][1]);  Bf[1][1] += reps;
  Bf[1][2] = myfabs(B[1][2]);  Bf[1][2] += reps;
  Bf[2][0] = myfabs(B[2][0]);  Bf[2][0] += reps;
  Bf[2][1] = myfabs(B[2][1]);  Bf[2][1] += reps;
  Bf[2][2] = myfabs(B[2][2]);  Bf[2][2] += reps;

  
#if TRACE1
  printf("Box test: Bf[3][3], B[3][3], T[3], a[3], b[3]\n");
  Mprintg(Bf);
  Mprintg(B);
  Vprintg(T);
  Vprintg(a);
  Vprintg(b);
#endif
  
  // if any of these tests are one-sided, then the polyhedra are disjoint
  r = 1;

  // A1 x A2 = A0
  t = myfabs(T[0]);
  
  r &= (t <= 
	  (a[0] + b[0] * Bf[0][0] + b[1] * Bf[0][1] + b[2] * Bf[0][2]));
  if (!r) return 1;
  
  // B1 x B2 = B0
  s = T[0]*B[0][0] + T[1]*B[1][0] + T[2]*B[2][0];
  t = myfabs(s);

  r &= ( t <=
	  (b[0] + a[0] * Bf[0][0] + a[1] * Bf[1][0] + a[2] * Bf[2][0]));
  if (!r) return 2;
    
  // A2 x A0 = A1
  t = myfabs(T[1]);
  
  r &= ( t <= 
	  (a[1] + b[0] * Bf[1][0] + b[1] * Bf[1][1] + b[2] * Bf[1][2]));
  if (!r) return 3;

  // A0 x A1 = A2
  t = myfabs(T[2]);

  r &= ( t <= 
	  (a[2] + b[0] * Bf[2][0] + b[1] * Bf[2][1] + b[2] * Bf[2][2]));
  if (!r) return 4;

  // B2 x B0 = B1
  s = T[0]*B[0][1] + T[1]*B[1][1] + T[2]*B[2][1];
  t = myfabs(s);

  r &= ( t <=
	  (b[1] + a[0] * Bf[0][1] + a[1] * Bf[1][1] + a[2] * Bf[2][1]));
  if (!r) return 5;

  // B0 x B1 = B2
  s = T[0]*B[0][2] + T[1]*B[1][2] + T[2]*B[2][2];
  t = myfabs(s);

  r &= ( t <=
	  (b[2] + a[0] * Bf[0][2] + a[1] * Bf[1][2] + a[2] * Bf[2][2]));
  if (!r) return 6;

  // A0 x B0
  s = T[2] * B[1][0] - T[1] * B[2][0];
  t = myfabs(s);
  
  r &= ( t <= 
	(a[1] * Bf[2][0] + a[2] * Bf[1][0] +
	 b[1] * Bf[0][2] + b[2] * Bf[0][1]));
  if (!r) return 7;
  
  // A0 x B1
  s = T[2] * B[1][1] - T[1] * B[2][1];
  t = myfabs(s);

  r &= ( t <=
	(a[1] * Bf[2][1] + a[2] * Bf[1][1] +
	 b[0] * Bf[0][2] + b[2] * Bf[0][0]));
  if (!r) return 8;

  // A0 x B2
  s = T[2] * B[1][2] - T[1] * B[2][2];
  t = myfabs(s);

  r &= ( t <=
	  (a[1] * Bf[2][2] + a[2] * Bf[1][2] +
	   b[0] * Bf[0][1] + b[1] * Bf[0][0]));
  if (!r) return 9;

  // A1 x B0
  s = T[0] * B[2][0] - T[2] * B[0][0];
  t = myfabs(s);

  r &= ( t <=
	  (a[0] * Bf[2][0] + a[2] * Bf[0][0] +
	   b[1] * Bf[1][2] + b[2] * Bf[1][1]));
  if (!r) return 10;

  // A1 x B1
  s = T[0] * B[2][1] - T[2] * B[0][1];
  t = myfabs(s);

  r &= ( t <=
	  (a[0] * Bf[2][1] + a[2] * Bf[0][1] +
	   b[0] * Bf[1][2] + b[2] * Bf[1][0]));
  if (!r) return 11;

  // A1 x B2
  s = T[0] * B[2][2] - T[2] * B[0][2];
  t = myfabs(s);

  r &= (t <=
	  (a[0] * Bf[2][2] + a[2] * Bf[0][2] +
	   b[0] * Bf[1][1] + b[1] * Bf[1][0]));
  if (!r) return 12;

  // A2 x B0
  s = T[1] * B[0][0] - T[0] * B[1][0];
  t = myfabs(s);

  r &= (t <=
	  (a[0] * Bf[1][0] + a[1] * Bf[0][0] +
	   b[1] * Bf[2][2] + b[2] * Bf[2][1]));
  if (!r) return 13;

  // A2 x B1
  s = T[1] * B[0][1] - T[0] * B[1][1];
  t = myfabs(s);

  r &= ( t <=
	  (a[0] * Bf[1][1] + a[1] * Bf[0][1] +
	   b[0] * Bf[2][2] + b[2] * Bf[2][0]));
  if (!r) return 14;

  // A2 x B2
  s = T[1] * B[0][2] - T[0] * B[1][2];
  t = myfabs(s);

  r &= ( t <=
	  (a[0] * Bf[1][2] + a[1] * Bf[0][2] +
	   b[0] * Bf[2][1] + b[1] * Bf[2][0]));
  if (!r) return 15;

  return 0;  // should equal 0
}

#endif
