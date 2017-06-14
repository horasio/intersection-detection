/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */

#include "kmeans.h"
#include <algorithm>
#include <iostream>

using namespace std;

void kMeansInit(vector<Vec3f> & pts, const int start, const int end) {
  int a = start + ( lrand48() % ( end - start ) );
  Vec3f va = pts[a];
  sort(pts.begin()+start, pts.begin()+end, [&](const Vec3f & left, const Vec3f & right) {
        return (left-va).norm2() < (right-va).norm2();
      });
  float total(0.0f);
  for( int i = start+1; i < end; ++i )
    total += (pts[i]-va).norm2();
  float dice = drand48() * total;
  //float dice = (0.05f + drand48() * 0.95f) * total;
  total = 0.0f;
  int i;
  for( i = start+1; i < end; ++i ) {
    total += (pts[i]-va).norm2();
    if( total >= dice )
      break;
  }
  if( i == end ) cerr << "Holy cows!";
  if( start+1 < i )
    swap(pts[start+1], pts[i]);
  //cerr << "k-Means init for " << (end - start) << " vertice --> choose the " << (i-start) << "-th.";
}

void kMeans(vector<Vec3f> & pts, const int start, const int end, int & middle) {
  if( end - start < 2 ) {
    cerr << "BIG ERROR K-MEANS NOT ENOUGH POINTS";
    return;
  }
  if( end - start == 2 ) {
    middle = start + 1;
    return;
  }
  kMeansInit(pts, start, end); // stores initials seed in the first two positions.
  Vec3f seed[2] = { pts[start], pts[start+1] };
  bool changed(true);
  middle = start;
  int count(0);
  while( changed ) {
    ++count;
    int newMiddle(end);
    changed = false;
    for( int i = start; i < newMiddle; ) {
      float d0 = (pts[i] - seed[0]).norm2();
      float d1 = (pts[i] - seed[1]).norm2();
      if( d0 < d1 ) {
        ++i;
      } else {
        if( i < middle ) changed = true;
        swap(pts[i], pts[newMiddle-1]);
        --newMiddle;
      }
    }
    if( newMiddle != middle ) changed = true;
    middle = newMiddle;
    if( changed ) {
      if( middle > start ) {
        seed[0].set(0.0f, 0.0f, 0.0f);
        for( int i = start; i < middle; ++i )
          seed[0] = seed[0] + pts[i];
        seed[0] = (1.0f/(middle-start)) * seed[0];
      }
      if( end > middle ) {
        seed[1].set(0.0f, 0.0f, 0.0f);
        for( int i = middle; i < end; ++i )
          seed[1] = seed[1] + pts[i];
        seed[1] = (1.0f/(end - middle)) * seed[1];
      }
    }
  }
  if( middle == start ) cerr << "BIG ERROR K-MEANS A";
  if( end  ==  middle ) cerr << "BIG ERROR K-MEANS B";
  //std::cerr << " <" << count << " iterations> " << (middle-start) << " / " << (end-middle) << endl;
}
