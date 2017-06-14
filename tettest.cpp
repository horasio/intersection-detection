/* Author: Samuel Hornus <samuel.hornus@inria.fr>
 * Copyright © Inria, 2017
 * Licence: Creative Commons CC BY-ND 3.0 license available online at
 * http://creativecommons.org/licenses/by-nd/3.0/
 */

#include "Statistics.h"
#include "chronograph.h"
#include "spherical.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <sqlite3.h>

//#define WITH_CGAL
#ifdef WITH_CGAL
#include "CGALStuff.h"
#endif

#include "convexes.h"
#include "obb.h"
#include "gjk.h"

#define INTER_MAX_ITER 100

using namespace std;

typedef set<int> IntSet;

typedef Statistics<unsigned long> IntStat;

static TCLAP::CmdLine gCmdLine("Hi! I am Tettest and I like to intersect convexes! <3", ' ', "0.0");

static TCLAP::ValueArg<string> benchArg("b", "bench", "path to SQLite3 database file for dumping statistics", false, "", "path");

static TCLAP::SwitchArg tetArg("t", "tet", "Intersect pairs of tetrahedra", false);
static TCLAP::SwitchArg generalArg("g", "general", "Intersect pairs of convex polytopes", false);
static TCLAP::SwitchArg frustumSphereArg("", "frustum-sphere", "frustum/sphere intersections", false);
static TCLAP::SwitchArg frustumAABBArg("", "frustum-aabb", "frustum/aabb intersections", false);
static TCLAP::SwitchArg obbArg("o", "obb", "Intersect pairs of OBB", false);

static vector<string> tetMethods{"sat", "naive", "sphere", "sphered", "gjk", "gjkd", "hybrid", "greene", "cmpd"};
static TCLAP::ValuesConstraint<string> allowedAlgos(tetMethods);

static TCLAP::ValueArg<string> algoArg("a", "algorithm", "Intersection algorithm to use", false, "sphere", &allowedAlgos);

static TCLAP::SwitchArg sortedArg("s", "sorted-vertices", "Organize vertices into hierarchy", false);

static TCLAP::ValueArg<int> nboArg("n", "number-of-convexes", "Number of convexes", false, 1000, "positive integer");

static TCLAP::ValueArg<int> nbvArg("v", "number-of-vertices", "Number of vertices per convex", false, 10, "positive integer");

static TCLAP::ValueArg<float> spreadArg("p", "spread", "Maximum shift", false, 0.0, "positive real");

static int gNbV(0);
static long gNumHits(0);
static IntStat planeStat;
static float gSpread(0.0f);

// =========================================================

sqlite3 *db;
void openDatabase() { // ------------------ Database creation
  bool dumpStats = benchArg.getValue().length() > 0;
  if( ! dumpStats ) return;
  int dbret;
  dbret = sqlite3_open(benchArg.getValue().c_str(), &db);
  if( dbret ) {
    cerr << "\nCan't open database: %s" << sqlite3_errmsg(db);
    sqlite3_close(db); benchArg.reset(); return;
  }
  const char * SQLTableCreation = "CREATE TABLE IF NOT EXISTS STATISTICS("  \
                                   "DATE            TEXT    NOT NULL," \
                                   "N_OBJECT        INT     NOT NULL," \
                                   "N_PAIR          INT     NOT NULL," \
                                   "N_HITS          INT     NOT NULL," \
                                   "N_VERTEX        INT     NOT NULL," \
                                   "ALGORITHM       TEXT    NOT NULL," \
                                   "SORTED          INT     NOT NULL," \
                                   "TEST_NAME       TEXT    NOT NULL," \
                                   "TIME            REAL    NOT NULL," \
                                   "SPREAD          INT     NOT NULL," \
                                   "N_FAIL          INT     NOT NULL," \
                                   "N_PLANE_MAX     INT     NOT NULL," \
                                   "N_PLANE         INT     NOT NULL," \
                                   "N_SQ_PLANE      INT     NOT NULL);";
  sqlite3_stmt * statement;
  const char * pzTail;
  dbret = sqlite3_prepare_v2(db, SQLTableCreation, strlen(SQLTableCreation), &statement, &pzTail);
  if( NULL == statement ) {
    cerr << "\nCan't prepare statement for table creation: " << sqlite3_errmsg(db)
      << ". retcode = " << dbret;
    sqlite3_close(db); benchArg.reset(); return;
  }
  dbret = sqlite3_step(statement);
  if( dbret != SQLITE_DONE ) {
    cerr << "\nCan't execute table creation: %s" << sqlite3_errmsg(db);
    sqlite3_close(db); benchArg.reset(); return;
  }
  sqlite3_finalize(statement);
  statement = NULL;
  return;
}

Chronograph myChrono;
unsigned long planeStatPerPair;
unsigned long nbFails;

void closeDatabase() {
  bool dumpStats = benchArg.getValue().length() > 0;
  if( ! dumpStats ) return;
  stringstream s;
  s.precision(12);
  s << "insert into statistics values (datetime('now','localtime'),"
    << nboArg.getValue() << ','
    << planeStat.nbSamples() << ','
    << gNumHits << ','
    << gNbV << ','
    << '\'' << (algoArg.getValue()) << "',"
    << (sortedArg.getValue() ? 1 : 0) << ','
    << ( tetArg.getValue() ? "'TETS'" :
        ( obbArg.getValue() ? "'OBB'" :
          ( generalArg.getValue() ? "'GENERAL'" :
            ( frustumSphereArg.getValue() ? "'FRUSTUM-SPHERE'" :
              ( frustumAABBArg.getValue() ? "'FRUSTUM-AABB'" : "'UNKNOWN'"
              )
            )
          )
        )
       ) << ','
    << myChrono.elapsed_time() << ','
    << static_cast<int>(gSpread*100.0f) << ','
    << nbFails << ','
    << (planeStat.max() == INTER_MAX_ITER ? planeStat.second_max() : planeStat.max()) << ','
    << planeStat.total() << ','
    << planeStat.sqTotal()
    << ");";
  char * err(NULL);
  sqlite3_exec(db, s.str().c_str(), NULL, NULL, &err);
  if( NULL != err )
    cerr << endl << err;
  sqlite3_close(db);
}

// =========================================================

bool readCommandLine(int argc, char **argv) {
  try
  {
    vector<TCLAP::Arg*> v{ & tetArg, & generalArg, & frustumSphereArg, & frustumAABBArg, & obbArg };
    gCmdLine.xorAdd(v);
    gCmdLine.add(benchArg);
    gCmdLine.add(algoArg);
    gCmdLine.add(sortedArg);
    gCmdLine.add(nbvArg);
    gCmdLine.add(nboArg);
    gCmdLine.add(spreadArg);
    gCmdLine.parse(argc, argv);

    if( spreadArg.getValue() < 0.0f )
      throw TCLAP::ArgException("bad value", "-p --spread");
    gSpread = spreadArg.getValue();

    if( nboArg.getValue() < 2 || nboArg.getValue() > 10000 )
      throw TCLAP::ArgException("bad value", "-n --number-of-convexes");

    if( nbvArg.getValue() < 3 || nbvArg.getValue() > 10000 )
      throw TCLAP::ArgException("bad value", "-v --number-of-vertices");
    gNbV = nbvArg.getValue();

    if( generalArg.getValue() && algoArg.getValue() == "naive" && nbvArg.getValue() > 1000 )
      throw TCLAP::ArgException("bad value", "too many vertices. Will segfault by filling the stack.");

    return false;
  }
  catch (const TCLAP::ArgException & e)
  {
    cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
  }
  catch (...)  // catch any exceptions
  {
    cerr << "Error: unknown exception caught" << endl;
  }
  return true;
}

// =========================================================

static const int edgesIdx[6][2]   = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
static const int normalsIdx[6][2] = {{2,3},{1,3},{1,2},{0,3},{0,2},{0,1}};
/*
 * 00 11 22 33 44 55
 * 01 12 23 34 45 50
 * 02 13 24 35 40 51
 * 03 14 25 30 41 52
 * 04 15 20 31 42 53
 * 05 10 21 32 43 54
 * */

bool faceSeparation(const Tet & t1, const Tet & t2) {
  return t1.faceSeparate(t2) || t2.faceSeparate(t1);
}

// ----------------------------------------------------------- AABB / Frustum

bool gjkDisjointAAB(const OBB & a, const VerticesOnly & b) {
#define TT_GJK_MACRO \
  a.AAmaximizeInDirection(dir, maxiA, vA); \
  b.minimizeInDirection(dir, miniB, vB); \
  ++planeStatPerPair; \
  if( maxiA <= miniB ) return true
  static GJK simplex;
  int vB; Vec3f vA;
  float maxiA, miniB;
  Vec3f dir(b.vertex(0) - a.center());
  TT_GJK_MACRO;
  simplex.set(dir);
  do {
    if( simplex.add_and_nearest(b.vertex(vB) - vA, dir) ) return false;
    TT_GJK_MACRO;
    if(planeStatPerPair >= INTER_MAX_ITER) {
      //cerr << "MAX NUMBER OF ITERATIONS EXCEEDED\n";
      ++nbFails;
      return false;
    }
  } while( true );
#undef TT_GJK_MACRO
}

bool hybridDisjointAAB(const OBB & a, const VerticesOnly & b) {
#define TT_MACRO \
  a.AAmaximizeInDirection(dir, maxiA, vA); \
  b.minimizeInDirection(dir, miniB, vB); \
  ++planeStatPerPair; \
  if( maxiA <= miniB ) return true
  static GJK simplex;
  static SphericalPolygon positiveBound, tempPoly;
  int vB; Vec3f vA;
  float maxiA, miniB;
  Vec3f dir(b.vertex(0) - a.center());
  TT_MACRO;
  simplex.set(dir);
  for( int n = 1; n <= 3; ++n ) {
    if( simplex.add_and_nearest(b.vertex(vB) - vA, dir) ) return false;
    TT_MACRO;
  }
  positiveBound.clear();
  positiveBound.emplace_back(simplex.pts_[0]);
  positiveBound.clip(simplex.pts_[1], tempPoly); positiveBound.swap(tempPoly);
  if( simplex.n_ > 2 ) {
    positiveBound.clip(simplex.pts_[2], tempPoly); positiveBound.swap(tempPoly);
  }
  positiveBound.clip(b.vertex(vB) - vA, tempPoly); positiveBound.swap(tempPoly);
  /**/
  if( positiveBound.empty() ) return false;
  do {
    dir = positiveBound.averageDirection();
    TT_MACRO;
    if(planeStatPerPair >= INTER_MAX_ITER) {
      ++nbFails;
      return false;
    }
    positiveBound.clip(b.vertex(vB) - vA, tempPoly); positiveBound.swap(tempPoly);
    if( positiveBound.empty() ) return false;
  } while( true );
#undef TT_MACRO
}

bool sphericalDisjointAAB(const OBB & a, const VerticesOnly & b) {
#define TT_SPH_MACRO \
  a.AAmaximizeInDirection(dir, maxiA, vA); \
  b.minimizeInDirection(dir, miniB, vB); \
  ++planeStatPerPair; \
  if( maxiA <= miniB ) return true
  static SphericalPolygon positiveBound, tempPoly;
  Vec3f vA; int vB;
  float maxiA, miniB;
  Vec3f dir(b.vertex(0) - a.center());
  TT_SPH_MACRO;
  positiveBound.clear();
  positiveBound.emplace_back(dir);
  positiveBound.clip((b.vertex(vB) - vA), tempPoly); positiveBound.swap(tempPoly);
  if( positiveBound.empty() ) return false;
  do {
    dir = positiveBound.averageDirection();
    TT_SPH_MACRO;
    if(planeStatPerPair >= INTER_MAX_ITER) {
      ++nbFails;
      return false;
    }
    positiveBound.clip((b.vertex(vB) - vA), tempPoly); positiveBound.swap(tempPoly);
    if( positiveBound.empty() ) return false;
  } while( true );
#undef TT_SPH_MACRO
}

// -------------------------------------------------------- Centered AABB / OBB

bool gjkDisjointOBB(const OBB & a, const OBB & b) {
#define TT_GJK_MACRO(testdir) \
  a.centeredAAmaximizeInDirection(testdir, maxiA, vA); \
  b.minimizeInDirection(testdir, miniB, vB); \
  ++planeStatPerPair; \
  if( maxiA <= miniB ) return true
  static GJK simplex;
  Vec3f vA, vB;
  float maxiA, miniB;
  Vec3f dir(b.center());// - a.center()); because a.center() == (0,0,0)
  simplex.clear();
  TT_GJK_MACRO(dir);
  simplex.set(dir);
  do {
    if( simplex.add_and_nearest(vB - vA, dir) ) return false;
    TT_GJK_MACRO(dir);
    if(planeStatPerPair >= INTER_MAX_ITER) {
      ++nbFails;
      return false;
    }
  } while( true );
#undef TT_GJK_MACRO
}

bool sphericalDisjointOBB(const OBB & a, const OBB & b) {
#define TT_SPH_MACRO(testdir) \
  a.centeredAAmaximizeInDirection(testdir, maxiA, vA); \
  b.minimizeInDirection(testdir, miniB, vB); \
  ++planeStatPerPair; \
  if( maxiA <= miniB ) return true
  static SphericalPolygon positiveBound, tempPoly;
  Vec3f vA; float maxiA;
  Vec3f vB; float miniB;
  Vec3f dir;
  TT_SPH_MACRO(b.center());
  positiveBound.clear();
  positiveBound.emplace_back(b.center());
  positiveBound.clip(vB - vA, tempPoly); positiveBound.swap(tempPoly);
  if( positiveBound.empty() ) return false;
  do {
    dir = positiveBound.averageDirection();
    TT_SPH_MACRO(dir);
    if(planeStatPerPair >= INTER_MAX_ITER) {
      ++nbFails;
      return false;
    }
    positiveBound.clip(vB - vA, tempPoly); positiveBound.swap(tempPoly);
    if( positiveBound.empty() ) return false;
  } while( true );
#undef TT_SPH_MACRO
}

// ------------------------------------------------------------ Convex / Sphere

template< typename Convex >
bool gjkDisjointWithSphere(const Convex & a, const Sphere & s) {
#define TT_GJK_MACRO(testdir) \
  a.maximizeInDirection(testdir, maxiA, vA); \
  s.minimizeInDirection(testdir, miniS, vS); \
  ++planeStatPerPair; \
  if( maxiA <= miniS ) return true
  static GJK simplex;
  int vA; Vec3f vS;
  float maxiA, miniS;
  Vec3f dir(s.vertex(0) - a.vertex(0));
  simplex.clear();
  TT_GJK_MACRO(dir);
  simplex.set(dir);
  do {
    if( simplex.add_and_nearest(vS - a.vertex(vA), dir) ) return false;
    TT_GJK_MACRO(dir);
    if(planeStatPerPair >= INTER_MAX_ITER) {
      ++nbFails;
      return false;
    }
  } while( true );
#undef TT_GJK_MACRO
}

template< typename Convex >
bool sphericalDisjointWithSphere(const Convex & a, const Sphere & s) {
#define TT_SPH_MACRO(testdir) \
  a.maximizeInDirection(testdir, maxiA, vA); \
  s.minimizeInDirection(testdir, miniS, vS); \
  ++planeStatPerPair; \
  if( maxiA <= miniS ) return true
  static SphericalPolygon positiveBound, tempPoly;
  int vA; Vec3f vS;
  float maxiA, miniS;
  Vec3f dir(s.vertex(0) - a.vertex(0));
  TT_SPH_MACRO(dir);
  positiveBound.clear();
  positiveBound.emplace_back(dir);
  positiveBound.clip(vS - a.vertex(vA), tempPoly); positiveBound.swap(tempPoly);
  if( positiveBound.empty() ) return false;
  do {
    dir = positiveBound.averageDirection();
    TT_SPH_MACRO(dir);
    if(planeStatPerPair >= INTER_MAX_ITER) {
      ++nbFails;
      return false;
    }
    positiveBound.clip(vS - a.vertex(vA), tempPoly); positiveBound.swap(tempPoly);
    if( positiveBound.empty() ) return false;
  } while( true );
#undef TT_SPH_MACRO
}

bool hybridDisjointWithSphere(const VerticesOnly & f, const Sphere & s) {
#define TT_MACRO \
  f.maximizeInDirection(dir, maxiF, vF); \
  s.minimizeInDirection(dir, miniS, vS); \
  ++planeStatPerPair; \
  if( maxiF <= miniS ) return true
  static GJK simplex;
  static SphericalPolygon positiveBound, tempPoly;
  int vF; Vec3f vS;
  float maxiF, miniS;
  Vec3f dir(s.center() - f.vertex(0));
  TT_MACRO;
  simplex.set(dir);
  for( int n = 1; n <= 3; ++n ) {
    if( simplex.add_and_nearest(vS - f.vertex(vF), dir) ) return false;
    TT_MACRO;
  }
  positiveBound.clear();
  positiveBound.emplace_back(simplex.pts_[0]);
  positiveBound.clip(simplex.pts_[1], tempPoly); positiveBound.swap(tempPoly);
  if( simplex.n_ > 2 ) {
    positiveBound.clip(simplex.pts_[2], tempPoly); positiveBound.swap(tempPoly);
  }
  positiveBound.clip(vS - f.vertex(vF), tempPoly); positiveBound.swap(tempPoly);
  /**/
  if( positiveBound.empty() ) return false;
  do {
    dir = positiveBound.averageDirection();
    TT_MACRO;
    if(planeStatPerPair >= INTER_MAX_ITER) {
      ++nbFails;
      return false;
    }
    positiveBound.clip(vS - f.vertex(vF), tempPoly); positiveBound.swap(tempPoly);
    if( positiveBound.empty() ) return false;
  } while( true );
#undef TT_MACRO
}

// -------------------------------------------------------------------- Convex / Convex

template< class Convex >
bool gjkDisjoint(const Convex & a, const Convex & b) {
  static GJK simplex;
  int vA, vB;
  Vec3f dir = b.vertex(0) - a.vertex(0);
  simplex.clear();
  if( ! a.differenceCoversZeroInDir(b, vA, vB, dir) ) return true;
  simplex.set(dir);
  do {
    if( simplex.add_and_nearest(b.vertex(vB) - a.vertex(vA), dir) ) return false;
    if( ! a.differenceCoversZeroInDir(b, vA, vB, dir) ) return true;
    if(planeStatPerPair >= INTER_MAX_ITER) {
      ++nbFails;
      return false;
    }
  } while( true );
}

template< class Convex >
float gjkDistance(const Convex & a, const Convex & b) {
  static GJK simplex;
  int vA(0), vB(0);
  Vec3f dir;
  simplex.clear();
  float upperBound(std::numeric_limits<float>::max());
  float lowerBound(0.0f);
  do {
    float maxOverA, minOverB;
    if( simplex.add_and_nearest_for_distance(b.vertex(vB) - a.vertex(vA), dir) ) return 0.0f;
    float d = dir.length();
    upperBound = min(upperBound, d);
    dir = dir / d;
    a.differenceInDir(b, vA, vB, maxOverA, minOverB, dir);
    lowerBound = max(lowerBound, (minOverB - maxOverA));
    if( upperBound - lowerBound < - 1e-5f ) {
      ++nbFails;
      cerr << "BUG: ";
    }
    if( upperBound - lowerBound < 1e-4f ) {
      return 0.5f * (upperBound + lowerBound);
    }
    if(planeStatPerPair >= INTER_MAX_ITER) {
      ++nbFails;
      return 0.5f * (upperBound + lowerBound);
    }
  } while( true );
}

//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

static SphericalPolygon tupTemp_;
struct ThreeUniquePoints {
  vector<Vec3f> pts_;
  int pos_ = 0;
  ThreeUniquePoints() : pts_() {}
  void add(const Vec3f & v) {
    bool same_found(false);
    for( const auto & p : pts_ ) {
      if( v == p ) {
        same_found = true;
        break;
      }
    }
    if( same_found ) return;
#define TUPSIZE 3
    if( pts_.size() < TUPSIZE ) {
      pts_.push_back(v);
      return;
    } else {
      pts_[pos_] = v;
      pos_ = (pos_ + 1) % TUPSIZE;
    }
  }
  void clip(const SphericalPolygon & src, SphericalPolygon & dest, const Vec3f & delta) {
    int i(1);
    src.clip(pts_[0] - delta, dest);
    for( ; i < pts_.size(); ) {
      dest.clip(pts_[i++] - delta, tupTemp_);
      dest.swap(tupTemp_);
    }
  }
};

//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

template< typename Convex >
bool sphericalDisjoint(const Convex & a, const Convex & b) {
  static SphericalPolygon positiveBound, tempPoly;
  int vA, vB;
  Vec3f dir(b.vertex(0) - a.vertex(0));
  if( ! a.differenceCoversZeroInDir(b, vA, vB, dir) ) return true;
  positiveBound.clear();
  positiveBound.emplace_back(dir);
  positiveBound.clip(b.vertex(vB) - a.vertex(vA), tempPoly); positiveBound.swap(tempPoly);
  if( positiveBound.empty() ) return false;
  do {
    if( ! a.differenceCoversZeroInDir(b, vA, vB, positiveBound.averageDirection()) ) return true;
    if(planeStatPerPair >= INTER_MAX_ITER) {
      ++nbFails;
      return false;
    }
    positiveBound.clip(b.vertex(vB) - a.vertex(vA), tempPoly); positiveBound.swap(tempPoly);
    if( positiveBound.empty() ) return false;
  } while( true );
}

//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

template< typename Convex >
bool sphericalDisjointWithBounds(const Convex & a, const Convex & b,
    SphericalPolygon & positiveBound, GJK & simplex,
    float & lowerBound, float & upperBound) {
    //Vec3f & P1, Vec3f & P2) {
  static SphericalPolygon tempPoly;
  static vector<Vec3f> vertices;
  int vA, vB;
  Vec3f P(b.vertex(0) - a.vertex(0));
  simplex.clear();
  upperBound = simplex.add_and_distance(P);
  positiveBound.clear();
  positiveBound.emplace_back(P);
  do {
    float maxOverA, minOverB;
    Vec3f n = positiveBound.averageDirection().normalized();
    a.differenceInDir(b, vA, vB, maxOverA, minOverB, n);
    if( maxOverA < minOverB ) {
      P = b.vertex(vB) - a.vertex(vA);
      upperBound = std::min(upperBound, simplex.add_and_distance(P));
      //upperBound = std::min(upperBound, P.norm2());
      lowerBound = minOverB - maxOverA; // > 0

      Vec3f delta(lowerBound * n);
      cout << endl << "NO INTERSECTION. low = " << lowerBound << ", hi = " << upperBound << endl;
      positiveBound.clip(P, P - delta, tempPoly); positiveBound.swap(tempPoly);
      /*vertices.clear();
      for( const auto & v : positiveBound )
        vertices.push_back(v.silVertex_);
      for( const auto & v : vertices ) {
        positiveBound.clip(v, v - delta, tempPoly, false); positiveBound.swap(tempPoly);
        if( positiveBound.empty() ) break;
      }*/
      return true;
    }
    if(planeStatPerPair >= INTER_MAX_ITER) {
      ++nbFails;
      return false;
    }
    P = b.vertex(vB) - a.vertex(vA);
    positiveBound.clip(P, tempPoly); positiveBound.swap(tempPoly);
    if( positiveBound.empty() ) return false;
    upperBound = std::min(upperBound, simplex.add_and_distance(P));
    //upperBound = std::min(upperBound, P.norm2());
  } while( true );
}

//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
//
//  TODO : reput ThreeUniquePoints (lasts) but use it to compute distance to
//  the triangle formed by the last 3 points, in order to have an upper bound

template< typename ConvexA, typename ConvexB >
float sphericalDistance(const ConvexA & a, const ConvexB & b) {
  static SphericalPolygon S_low, S_cur, tempPoly;
  static GJK simplex;
  float low, low_diff;
  Vec3f low_n, low_p;
  float cur, cur_diff;
  Vec3f cur_n, cur_p, dir;
  float hi(std::numeric_limits<float>::max());
  float maxOverA, minOverB;
  int vA, vB;
  if( ! sphericalDisjointWithBounds(a, b, S_cur, simplex, low, hi) )
    return 0.0f;
  //hi = ::sqrtf(hi);
#define TRACE
  if( S_cur.empty() ) {
#ifdef TRACE
    cout << endl << "EMPTY. " << (hi-low) << ", low:"<<low<<", cur:"<<cur<<", hi:"<<hi << ", |S|:" << S_cur.size() << endl;
#endif
    return low;
  }
  if( hi - low < 1e-4f ) {
#ifdef TRACE
    cout << endl << "GOOD BOUNDS. " << (hi-low) << ", low:"<<low<<", cur:"<<cur<<", hi:"<<hi << ", |S|:" << S_cur.size() << endl;
#endif
    return 0.5f*(low+hi);
  }

  S_low = S_cur;

  float k = 2.0f;
  float kk = 2.0f*k-1.0f;

  low_n = S_cur.averageDirection().normalized();
  a.differenceInDir(b, vA, vB, maxOverA, minOverB, low_n);
  low_diff = minOverB - maxOverA;
  low_p = b.vertex(vB) - a.vertex(vA);
  hi = std::min(hi, simplex.add_and_distance(low_p));

  cur_n = low_n;
  cur_diff = low_diff;
  cur_p = low_p;

  cur = (kk * low < hi) ? k * low : 0.5f * (low + hi);

  bool cur_is_low(true);

  do {
#ifdef TRACE
    cout << endl << (hi-low) << ", low:"<<low<<", cur:"<<cur<<", hi:"<<hi << ", |S|:" << S_cur.size();
#endif
    if( hi - low < 1e-4f ) {
#ifdef TRACE
      cout << " FOUND! " << 0.5f*(low+hi)<<". ";
#endif
      return 0.5f*(low+hi);
    }
    if( planeStatPerPair >= INTER_MAX_ITER ) {
      ++nbFails;
#ifdef TRACE
      cout << " FAIL " << 0.5f*(low+hi)<<". ";
#endif
      return 0.5f*(low+hi);
    }

    if( cur_diff <= cur ) { // maxOverA + cur >= minOverB
      S_cur.clip(cur_p, cur_p - cur * cur_n, tempPoly);// S_cur.swap(tempPoly);
      if( tempPoly.empty() ) { // cur is too large, restart on S_low
        hi = cur;
        cur = (kk * low < hi) ? k * low : 0.5f * (low + hi);
        if( ! cur_is_low ) {
          S_cur = S_low;
          cur_is_low = true;
          cur_diff = low_diff;
          cur_n = low_n;
          cur_p = low_p;
        }
#ifdef TRACE
        cout<<" --> too large, restart on S_low. ";
#endif
      } else { // continue clipping
        S_cur.swap(tempPoly);
        cur_is_low = false;
        cur_n = S_cur.averageDirection().normalized();
        a.differenceInDir(b, vA, vB, maxOverA, minOverB, cur_n);
        cur_diff = minOverB - maxOverA;
        cur_p = b.vertex(vB) - a.vertex(vA);
        hi = std::min(hi, simplex.add_and_distance(cur_p));
        if( hi <= cur ) {
          cur = (kk * low < hi) ? k * low : 0.5f * (low + hi);
          S_cur = S_low;
          cur_is_low = true;
          cur_diff = low_diff;
          cur_n = low_n;
          cur_p = low_p;
        }
#ifdef TRACE
        cout << " --> ...";
#endif
      }
    } else { // maxOverA + cur < minOverB => new lower bound: cur_diff > cur
      Vec3f delta(cur_diff * cur_n);
      S_cur.clip(cur_p, cur_p - delta, S_low);
      /*for( int i = 0; i < S_cur.size(); ++i ) {
        tempPoly.swap(S_low);
        tempPoly.clip(tempPoly[i].silVertex_, tempPoly[i].silVertex_ - delta, S_low, false);
        if( S_low.empty() ) break;
      }*/
      if( S_low.empty() ) {
#ifdef TRACE
        cout << " not enough precision: "<<cur_diff <<". ";
#endif
        return cur_diff;
      }
      //S_cur.swap(tempPoly);
      //S_low = S_cur;
      S_cur = S_low;
      cur_is_low = true;
      low = cur_diff;
      low_n = S_cur.averageDirection().normalized();
      a.differenceInDir(b, vA, vB, maxOverA, minOverB, low_n);
      low_diff = minOverB - maxOverA;
      low_p = b.vertex(vB) - a.vertex(vA);
      hi = std::min(hi, simplex.add_and_distance(low_p));
      cur_diff = low_diff;
      cur_n = low_n;
      cur_p = low_p;
      //hi = std::min(hi, low_p.length());
      cur = (kk * low < hi) ? k * low : 0.5f * (low + hi);
#ifdef TRACE
      cout << " --> too low, new lower bound";
#endif
    }

  } while( true );
}

// =============================================================

Vec3f genConvex(int nbVertices, Polyhedron_3 & poly) {
  if( nbVertices < 4 )
    nbVertices = 4;
  static vector<Point_3> v;
  v.clear();
  RandomPointOnSphere rps;
  Vec3f avg;
  float shift = drand48() * gSpread;
  Vec3f s(shift, 0.0f, 0.0f);
  for( int i = 0; i < nbVertices; ++i ) {
    Vec3f p(rps.vec3f() + s);
    v.emplace_back(p.x(), p.y(), p.z());
    avg = avg + p;
  }
  avg = (1.0f/nbVertices) * avg;
  convex_hull_3(v.begin(), v.end(), poly);
  return avg;
}

// =============================================================

template< typename Convex >
bool testDisjoint(const Convex & a, const Convex & b) {
  // test B against facets of A
  const int nva = a.nbVertices();
  const int nvb = b.nbVertices();
  const int nfa = a.nbFacets();
  const int nfb = b.nbFacets();
  for( int af = 0; af < nfa; ++af ) {
    const Vec3f & n = a.facets_[af].normal_;
    float minOverB(n | b.vertices_[0]);
    for( int bv = 1; bv < nvb; ++bv ) {
      float d = b.vertices_[bv] | n;
      if( minOverB > d ) minOverB = d;
      if( a.facets_[af].maxAlong_ > minOverB ) break;
    }
    ++planeStatPerPair;
    if( a.facets_[af].maxAlong_ <= minOverB ) return true;
  }
  float AvDotBn[nva][nfb];
  // test A against facets of B
  for( int bf = 0; bf < nfb; ++bf ) {
    const Vec3f & n = b.facets_[bf].normal_;
    float minOverA(n | a.vertices_[0]);
    AvDotBn[0][bf] = minOverA;
    for( int av = 1; av < nva; ++av ) {
      float d;
      AvDotBn[av][bf] = d = a.vertices_[av] | n;
      if( minOverA > d ) minOverA = d;
    }
    ++planeStatPerPair;
    if( b.facets_[bf].maxAlong_ <= minOverA ) return true;
  }
  // test A-B-edge pairs
  float bpos[nfb];
  const int nea = a.nbEdges();
  const int neb = b.nbEdges();
  for( int ia = 0; ia < nea; ++ia ) {
    int ani0 = a.edges_[ia][2];
    int ani1 = a.edges_[ia][3];
    for( int i = 0; i < nfb; ++i ) {// precompute the dot-product of normals of |b| with |a|'s current edge
      bpos[i] = AvDotBn[a.edges_[ia][0]][i] - AvDotBn[a.edges_[ia][1]][i];
    }
    for( int ib = 0; ib < neb; ++ib ) {
      int bni0 = b.edges_[ib][2];
      int bni1 = b.edges_[ib][3];
      if( bpos[bni0] * bpos[bni1] >= 0.0f ) continue;
      float bpos0 = ::fabsf(bpos[bni0]);
      float bpos1 = ::fabsf(bpos[bni1]);
      float maxOverB = bpos0 * b.facets_[bni1].maxAlong_ + bpos1 * b.facets_[bni0].maxAlong_;
      float xbDotAw0 = bpos0 * AvDotBn[ani0][bni1] + bpos1 * AvDotBn[ani0][bni0];
      float xbDotAw1 = bpos0 * AvDotBn[ani1][bni1] + bpos1 * AvDotBn[ani1][bni0];
      int edId = a.edges_[ia][0];
      if( xbDotAw0 > xbDotAw1 ) xbDotAw0 = xbDotAw1;
      float xbDotAve = bpos0 * AvDotBn[edId][bni1] + bpos1 * AvDotBn[edId][bni0];
      ++planeStatPerPair;
      if( maxOverB <= xbDotAve  && xbDotAve <= xbDotAw0 ) return true;
    }
  }
  return false;
}

template< >
bool testDisjoint<VerticesOnly>(const VerticesOnly & a, const VerticesOnly & b) {
  return true;
}

template< >
bool testDisjoint<SortedVertices>(const SortedVertices & a, const SortedVertices & b) {
  return true;
}

// =============================================================

template< typename Convex >
void test_general(const int N, const int nbv) {

  vector<Convex> convexes;
  convexes.reserve(N);

  int algo = 0;
  if( algoArg.getValue() == "naive" ) {
    algo = 0;
  } else if( algoArg.getValue() == "sphere" ) {
    algo = 1;
  } else if( algoArg.getValue() == "sat" ) {
    algo = 1;
    cerr << "WARNING: SAT not implemented for general convexes. Using 'sphere'";
  } else if( algoArg.getValue() == "gjk" ) {
    algo = 2;
  } else if( algoArg.getValue() == "sphered" ) {
    algo = 3;
  } else if( algoArg.getValue() == "gjkd" ) {
    algo = 4;
  } else if( algoArg.getValue() == "cmpd" ) {
    algo = 5;
  }
  if( algo == 0 ) {
    Polyhedron_3 poly;
    for( int i = 0; i < N; ++i ) {
      Vec3f center = genConvex(nbv, poly);
      convexes.emplace_back(poly);
      convexes.back().center_ = center;
    }
  } else {
    for( int i = 0; i < N; ++i ) {
      convexes.emplace_back(nbv, drand48() * gSpread);
    }
  }

  myChrono.start();

  for( int j = 1; j < N; ++j )
    for( int i = 0; i < j; ++i ) { // TEST GENERAL
      planeStatPerPair = 0;
      switch( algo ) {
        case 0: if( ! testDisjoint(convexes[i], convexes[j]) ) ++gNumHits; break;
        case 1: if( ! sphericalDisjoint(convexes[i], convexes[j]) ) ++gNumHits; break;
        case 2: if( ! gjkDisjoint(convexes[i], convexes[j]) ) ++gNumHits; break;
        case 3: if( 0.0f == sphericalDistance(convexes[i], convexes[j]) ) ++gNumHits; break;
        case 4: if( 0.0f == gjkDistance(convexes[i], convexes[j]) ) ++gNumHits; break;
        case 5: {
                  const float gjkd = gjkDistance(convexes[i], convexes[j]);
                  const float sphd = sphericalDistance(convexes[i], convexes[j]);
                  if( sphd > 0.0f && gjkd > 0.0f ) {
                    if( fabs(1.0f-sphd/gjkd) > 0.01f ) {
                      cerr << "** ";
                      cerr << "DSS=" << sphd << ". GJKD=" << gjkd << endl;
                    }
                  }
                  else
                    ++gNumHits;
                  break;
                }
        default: break;
      }
      planeStat.addSample(planeStatPerPair);
    }

  myChrono.stop();

  size_t pairs = planeStat.nbSamples();
  cout << gNumHits << " hits over " << pairs << " pairs (";
  cout << (100.0*(double)gNumHits)/pairs << "%) in " << myChrono.elapsed_time() << " seconds ( ";
  cout << (pairs / myChrono.elapsed_time() / 1000.0) << " KPairs/s ).";
  cout << endl << planeStat.min() << " / " << planeStat.mean() << "+/-" << planeStat.dev() << " / " << planeStat.second_max() << " / " << planeStat.max()
    << " min/avg/dev/premax/max directions/pair.";
  cout << endl << "FAILS:" << nbFails << endl << flush;
}

// =============================================================


int obbDisjoint(const OBB & a, const OBB & b) {
  Matrix3f R;
  Matrix::transMul(a.rotation_, b.rotation_, R);
  Vec3f T = a.rotation_.tMulVec3(b.center_ - a.center_);
  return obb_disjoint(R.d, T.data(), a.radii_.data(), b.radii_.data());
}

void test_obb(const int N) {
  if( gSpread == 0.0f )
    gSpread = 1e-2;
  vector<OBB> obbs;
  obbs.reserve(N);
  for( int j = 0; j < N; ++j ) {
    float x = drand48() * gSpread - gSpread / 2.0f;
    float y = drand48() * gSpread - gSpread / 2.0f;
    float z = drand48() * gSpread - gSpread / 2.0f;
    float dx = drand48();
    float dy = drand48();
    float dz = drand48();
    obbs.emplace_back(Vec3f(x,y,z), Vec3f(dx,dy,dz));
  }

  int algo = 0;
  if( algoArg.getValue() == "sat" )
    algo = 0;
  else if( algoArg.getValue() == "gjk" )
    algo = 2;
  else
    algo = 1;

  OBB o;
  myChrono.start();
  for( int j = 1; j < N; ++j ) {
    o.radii_ = obbs[j].radii_;
    for( int i = 0; i < j; ++i ) {
      planeStatPerPair = 0;
      switch( algo ) {
        case 0:
          planeStatPerPair = obbDisjoint(obbs[i], obbs[j]);
          if( planeStatPerPair == 0 ) {
            ++gNumHits; // in case of a hit, 15 SAT tests have been performed
            planeStatPerPair = 15;
          }
          break;
        case 1:
          Matrix::transMul(obbs[i].rotation_, obbs[j].rotation_, o.rotation_);
          o.center_ = obbs[i].rotation_.tMulVec3(obbs[j].center_ - obbs[i].center_);
          if( ! sphericalDisjointOBB(obbs[i], o) ) {
            ++gNumHits;
          }
          break;
        case 2:
          Matrix::transMul(obbs[i].rotation_, obbs[j].rotation_, o.rotation_);
          o.center_ = obbs[i].rotation_.tMulVec3(obbs[j].center_ - obbs[i].center_);
          if( ! gjkDisjointOBB(obbs[i], o) )
            ++gNumHits;
          break;
      }
      planeStat.addSample(planeStatPerPair);
    }
  }

  myChrono.stop();

  size_t pairs = planeStat.nbSamples();
  cout << gNumHits << " hits over " << pairs << " OBBs (";
  cout << (100.0*(double)gNumHits)/pairs << "%) in " << myChrono.elapsed_time() << " seconds ( ";
  cout << (pairs / myChrono.elapsed_time() / 1000.0) << " KPairs/s ).";
  cout << endl << planeStat.min() << " / " << planeStat.mean() << "+/-" << planeStat.dev() << " / " << planeStat.second_max() << " / " << planeStat.max()
    << " min/avg/dev/premax/max directions/pair.";
  cout << endl << flush;
}

// =============================================================

void test_aabb_frustum(const int N) {
  vector<OBB> boxes;
  boxes.reserve(N);
  Frustum f;
  float spr = f.radius() * (gSpread);
  for( int j = 0; j < N; ++j ) {
    float x, y, z;
    do {
      x = drand48() * 2.0f*spr - spr;
      y = drand48() * 2.0f*spr - spr;
      z = drand48() * 2.0f*spr - spr;
    } while( x*x+y*y+z*z > spr*spr);
    float dx = 0.0f+1.0f*drand48();
    float dy = 0.0f+1.0f*drand48();
    float dz = 0.0f+1.0f*drand48();
    boxes.emplace_back(Vec3f(x,y,z), Vec3f(dx,dy,dz));
  }
  if( algoArg.getValue() != "greene" ) {
    vector<VerticesOnly> frustums;
    for( int i = 0; i < N; ++i ) {
      frustums.emplace_back(1, 0.0f);
      frustums.back().makeFrustum();
    }
    int algo = 0; // sphere
    if( algoArg.getValue() == "gjk" )
      algo = 1;
    else if( algoArg.getValue() == "hybrid" )
      algo = 2;
    myChrono.start();

    for( int i = 0; i < N; ++i ) {
      for( int j = 0; j < N; ++j ) {
        planeStatPerPair = 0;
        switch( algo ) {
          case 0:
            if( ! sphericalDisjointAAB(boxes[i], frustums[j]) ) ++gNumHits; break;
          case 1:
            if( ! gjkDisjointAAB(boxes[i], frustums[j]) ) ++gNumHits; break;
          case 2:
            if( ! hybridDisjointAAB(boxes[i], frustums[j]) ) ++gNumHits; break;
          default: break;
        }
        planeStat.addSample(planeStatPerPair);
      }
    }
    myChrono.stop();
  } else { // Greene's algorithm
    vector<Frustum> frustums(N);
    myChrono.start();

    for( int i = 0; i < N; ++i ) {
      for( int j = 0; j < N; ++j ) {
        planeStatPerPair = 0;
        if( frustums[j].hitsAAB(boxes[i]) )
          ++gNumHits;
        planeStat.addSample(planeStatPerPair);
      }
    }
    myChrono.stop();
  }
  size_t pairs = planeStat.nbSamples();
  cout << gNumHits << " hits over " << pairs << " of frustum/aabb (";
  cout << (100.0*(double)gNumHits)/pairs << "%) in " << myChrono.elapsed_time() << " seconds ( ";
  cout << (pairs / myChrono.elapsed_time() / 1000.0) << " KPairs/s ).";
  cout << endl << planeStat.min() << " / " << planeStat.mean() << "+/-" << planeStat.dev() << " / " << planeStat.second_max() << " / " << planeStat.max()
    << " min/avg/dev/premax/max directions/pair.";
  cout << endl << flush;
}

// ============================= =================

void test_sphere_frustum(const int N) {
  vector<Sphere> spheres;
  spheres.reserve(N);
  Frustum f;
  float spr = f.radius() * (gSpread);
  for( int j = 0; j < N; ++j ) {
    float x, y, z;
    do {
      x = drand48() * 2.0f*spr - spr;
      y = drand48() * 2.0f*spr - spr;
      z = drand48() * 2.0f*spr - spr;
    } while( x*x+y*y+z*z > spr*spr);
    spheres.emplace_back(Vec3f(x,y,z), drand48());
  }
  vector<VerticesOnly> frustums;
  for( int i = 0; i < N; ++i ) {
    frustums.emplace_back(1, 0.0f);
    frustums.back().makeFrustum();
  }

  int algo = 0;
  if( algoArg.getValue() == "gjk" )
    algo = 1;
  else if( algoArg.getValue() == "hybrid" )
    algo = 2;

  myChrono.start();

  for( int j = 0; j < N; ++j )
    for( int i = 0; i < N; ++i ) {
      planeStatPerPair = 0;
      switch( algo ) {
        case 0:
          if( ! sphericalDisjointWithSphere(frustums[i], spheres[j]) ) ++gNumHits; break;
        case 1:
          if( ! gjkDisjointWithSphere(frustums[i], spheres[j]) ) ++gNumHits; break;
        case 2:
          if( ! hybridDisjointWithSphere(frustums[i], spheres[j]) ) ++gNumHits; break;
        default: break;
      }
      planeStat.addSample(planeStatPerPair);
    }

  myChrono.stop();

  size_t pairs = planeStat.nbSamples();
  cout << gNumHits << " hits over " << pairs << " of sphere/frustum (";
  cout << (100.0*(double)gNumHits)/pairs << "%) in " << myChrono.elapsed_time() << " seconds ( ";
  cout << (pairs / myChrono.elapsed_time() / 1000.0) << " KPairs/s ).";
  cout << endl << planeStat.min() << " / " << planeStat.mean() << "+/-" << planeStat.dev() << " / " << planeStat.second_max() << " / " << planeStat.max()
    << " min/avg/dev/premax/max directions/pair.";
  cout << endl << flush;
}

// ============================= =================

bool tetDisjointDirect(const Tet & a, const Tet & b) {
  // test B against facets of A
  for( int an = 0; an < 4; ++an ) {
    float minOverB(std::numeric_limits<float>::max());
    for( int bv = 0; bv < 4; ++bv ) {
      float d = b.vertices_[bv] | a.normals_[an];
      if( minOverB > d ) minOverB = d;
    }
    float maxOverA = a.vertices_[(an+1)&3] | a.normals_[an];
    if( maxOverA <= minOverB ) return true;
  }
  float AvDotBn[5][4];
  // test A against facets of B
  for( int bn = 0; bn < 4; ++bn ) {
    float minOverA(std::numeric_limits<float>::max());
    for( int av = 0; av < 4; ++av ) {
      float d;
      AvDotBn[av][bn] = d = a.vertices_[av] | b.normals_[bn];
      if( minOverA > d ) minOverA = d;
    }
    float maxOverB;
    AvDotBn[4][bn] = maxOverB = b.vertices_[(bn+1)&3] | b.normals_[bn];
    if( maxOverB <= minOverA ) return true;
  }
  // test A-B-edge pairs
  float bpos[4];
  for( int ia = 0; ia < 6; ++ia ) {
    int ani0 = normalsIdx[ia][0];
    int ani1 = normalsIdx[ia][1];
    for( int i = 0; i < 4; ++i ) {// precompute the dot-product of normals of |b| with |a|'s current edge
      bpos[i] = AvDotBn[edgesIdx[ia][0]][i] - AvDotBn[edgesIdx[ia][1]][i];
    }
    for( int ib = 0; ib < 6; ++ib ) {
      int bni0 = normalsIdx[ib][0];
      int bni1 = normalsIdx[ib][1];
      if( bpos[bni0] * bpos[bni1] >= 0.0f ) continue;
      float bpos0 = ::fabsf(bpos[bni0]);
      float bpos1 = ::fabsf(bpos[bni1]);
      float maxOverB = bpos0 * AvDotBn[ 4  ][bni1] + bpos1 * AvDotBn[ 4  ][bni0];
      float xbDotAw0 = bpos0 * AvDotBn[ani0][bni1] + bpos1 * AvDotBn[ani0][bni0];
      float xbDotAw1 = bpos0 * AvDotBn[ani1][bni1] + bpos1 * AvDotBn[ani1][bni0];
      int edId = edgesIdx[ia][0];
      if( xbDotAw0 > xbDotAw1 ) xbDotAw0 = xbDotAw1;
      float xbDotAve = bpos0 * AvDotBn[edId][bni1] + bpos1 * AvDotBn[edId][bni0];
      if( maxOverB <= xbDotAve  && xbDotAve <= xbDotAw0 ) return true;
    }
  }
  return false;
}

// --------------------------------------------------------

bool tetDisjoint(const Tet & a, const Tet & b) {
  // test B against facets of A
  for( int an = 0; an < 4; ++an ) {
    float minOverB(std::numeric_limits<float>::max());
    for( int bv = 0; bv < 4; ++bv ) {
      float d = b.vertices_[bv] | a.normals_[an];
      if( minOverB > d ) minOverB = d;
    }
    float maxOverA = a.vertices_[(an+1)&3] | a.normals_[an];
    ++planeStatPerPair;
    if( maxOverA <= minOverB ) return true;
  }
  float AvDotBn[5][4];
  // test A against facets of B
  for( int bn = 0; bn < 4; ++bn ) {
    float minOverA(std::numeric_limits<float>::max());
    for( int av = 0; av < 4; ++av ) {
      float d;
      AvDotBn[av][bn] = d = a.vertices_[av] | b.normals_[bn];
      if( minOverA > d ) minOverA = d;
    }
    float maxOverB;
    AvDotBn[4][bn] = maxOverB = b.vertices_[(bn+1)&3] | b.normals_[bn];
    ++planeStatPerPair;
    if( maxOverB <= minOverA ) return true;
  }
  // test A-B-edge pairs
  float bpos[4];
  for( int ia = 0; ia < 6; ++ia ) {
    for( int i = 0; i < 4; ++i ) {// precompute the dot-product of edges of |b| with |a|'s edge |ae|
      bpos[i] = AvDotBn[edgesIdx[ia][0]][i] - AvDotBn[edgesIdx[ia][1]][i];
    }
    for( int ib = 0; ib < 6; ++ib ) {
      int bni0 = normalsIdx[ib][0];
      int bni1 = normalsIdx[ib][1];
      if( bpos[bni0] * bpos[bni1] >= 0.0f ) continue;
      float bpos0 = ::fabsf(bpos[bni0]);
      float bpos1 = ::fabsf(bpos[bni1]);
      float maxOverBAlongBdir = bpos0 * AvDotBn[4][bni1] + bpos1 * AvDotBn[4][bni0];
      float v0, v1, v2, v3;
      v0 = bpos0 * AvDotBn[0][bni1] + bpos1 * AvDotBn[0][bni0];
      v1 = bpos0 * AvDotBn[1][bni1] + bpos1 * AvDotBn[1][bni0];
      if( v1 < v0 ) v0 = v1;
      v2 = bpos0 * AvDotBn[2][bni1] + bpos1 * AvDotBn[2][bni0];
      v3 = bpos0 * AvDotBn[3][bni1] + bpos1 * AvDotBn[3][bni0];
      if( v3 < v2 ) v2 = v3;
      if( v2 < v0 ) v0 = v2; 
      ++planeStatPerPair;
      if( v0 >= maxOverBAlongBdir )
        return true;
    }
  }
  return false;
}

bool edgeEdgeSeparationSAT(const Tet & t1, const Tet & t2) {
  Vec3f edges1[6], edges2[6];
  for( int i = 0; i < 6; ++i ) {
    edges1[i] = t1.vertices_[edgesIdx[i][0]]-t1.vertices_[edgesIdx[i][1]];
    edges2[i] = t2.vertices_[edgesIdx[i][0]]-t2.vertices_[edgesIdx[i][1]];
  }
  for( int e1 = 0; e1 < 6; ++e1 ) {
    for( int e2 = 0; e2 < 6; ++e2 ) {
      Vec3f direction = Vec3f::cross(edges1[e1], edges2[e2]);
      //if( direction.norm2() < 1e-6 ) continue;
      ++planeStatPerPair;
      if( t2.maxAlongDirection(direction) <= t1.minAlongDirection(direction) ) return true;
      if( t1.maxAlongDirection(direction) <= t2.minAlongDirection(direction) ) return true;
    }
  }
  return false;
}

bool tetDisjointSAT(const Tet & t1, const Tet & t2) {
  return faceSeparation(t1, t2) || edgeEdgeSeparationSAT(t1, t2);
}

#ifdef WITH_CGAL
bool tetDisjointCGAL(const Tet & t1, const Tet & t2) {
  Point_3 p10(t1.vertices_[0].x(), t1.vertices_[0].y(), t1.vertices_[0].z());
  Point_3 p11(t1.vertices_[1].x(), t1.vertices_[1].y(), t1.vertices_[1].z());
  Point_3 p12(t1.vertices_[2].x(), t1.vertices_[2].y(), t1.vertices_[2].z());
  Point_3 p13(t1.vertices_[3].x(), t1.vertices_[3].y(), t1.vertices_[3].z());
  Point_3 p20(t2.vertices_[0].x(), t2.vertices_[0].y(), t2.vertices_[0].z());
  Point_3 p21(t2.vertices_[1].x(), t2.vertices_[1].y(), t2.vertices_[1].z());
  Point_3 p22(t2.vertices_[2].x(), t2.vertices_[2].y(), t2.vertices_[2].z());
  Point_3 p23(t2.vertices_[3].x(), t2.vertices_[3].y(), t2.vertices_[3].z());
  Triangle_3 tri1[4]={{p10,p11,p12}, {p10,p11,p13}, {p10,p12,p13}, {p11,p12,p13}};
  Triangle_3 tri2[4]={{p20,p21,p22}, {p20,p21,p23}, {p20,p22,p23}, {p21,p22,p23}};
  Tetrahedron_3 tet1(p10,p11,p12,p13);
  Tetrahedron_3 tet2(p20,p21,p22,p23);
  for( int i = 0; i < 4; ++i )
    if( CGAL::do_intersect(tri1[i], tet2) )
      return false;
  for( int i = 0; i < 4; ++i )
    if( CGAL::do_intersect(tri2[i], tet1) )
      return false;
  return true;
}
#endif

bool positiveOrientation(const Tet & t) {
  Vec3f e1 = t.vertices_[1] - t.vertices_[0];
  Vec3f e2 = t.vertices_[2] - t.vertices_[0];
  Vec3f e3 = t.vertices_[3] - t.vertices_[0];
  return (Vec3f::cross(e1,e2) | e3) > 0.0f;
}

void genTet(Tet & t) {
recommence:
  float shift = drand48() * gSpread;
  RandomPointOnSphere rps;
  for( int i = 0; i < 4; ++i ) {
    t.vertices_[i] = rps.vec3f();
  }
  if( ! positiveOrientation(t) )
    swap(t.vertices_[0], t.vertices_[1]);
  t.computeNormalsAndCenter();
  if( ! t.containsZero() ) {
    goto recommence;
  }
  t.shift(shift);
#ifdef WITH_CGAL
  Point_3 p0(t.vertices_[0].x(), t.vertices_[0].y(), t.vertices_[0].z());
  Point_3 p1(t.vertices_[1].x(), t.vertices_[1].y(), t.vertices_[1].z());
  Point_3 p2(t.vertices_[2].x(), t.vertices_[2].y(), t.vertices_[2].z());
  Point_3 p3(t.vertices_[3].x(), t.vertices_[3].y(), t.vertices_[3].z());
  if( CGAL::orientation(p0, p1, p2, p3) != CGAL::POSITIVE )
    goto recommence;
#endif
}

size_t bothRight = 0;
size_t bothWrong = 0;
size_t satRight = 0;
size_t disRight = 0;

void printMistakes() {
  cout
    << endl << "           SAT right | SAT wrong"
    << endl << " SH right  " << bothRight << " | " << disRight
    << endl << " SH wrong  " << satRight << " | " << bothWrong
    << endl << flush;
}

int main(int argc, char **argv) {
  ::srand48(time(NULL));

  nbFails = 0;

  if( readCommandLine(argc, argv) )
    exit(-1);

  openDatabase();

  if( generalArg.getValue() ) {
    if( algoArg.getValue() == "naive" )
      test_general<FullConvex>(nboArg.getValue(), nbvArg.getValue());
    else if( sortedArg.getValue() )
      test_general<SortedVertices>(nboArg.getValue(), nbvArg.getValue());
    else
      test_general<VerticesOnly>(nboArg.getValue(), nbvArg.getValue());

  } else if( frustumAABBArg.getValue() ) {
    gNbV = 8;
    test_aabb_frustum(nboArg.getValue());

  } else if( frustumSphereArg.getValue() ) {
    gNbV = 8;
    test_sphere_frustum(nboArg.getValue());

  } else if( obbArg.getValue() ) {
    gNbV = 8;
    test_obb(nboArg.getValue());

  } else if( tetArg.getValue() ) { // TET TETRA
    gNbV = 4;
    const int N = nboArg.getValue();
    typedef vector<Tet> TetVector;
    TetVector tets(N);
    for( int i = 0; i < N; ++i ) {
      genTet(tets[i]);
    }
    long pairs(0);

    myChrono.start();

    int algo = 0;
    if( algoArg.getValue() == "sat" ) algo = 0;
    else if( algoArg.getValue() == "naive" ) algo = 1;
    else if( algoArg.getValue() == "sphere" ) algo = 2;
    else if( algoArg.getValue() == "gjk" ) algo = 3;
    else if( algoArg.getValue() == "sphered" ) algo = 4;
    else if( algoArg.getValue() == "gjkd" ) algo = 5;
    for( int j = N-1; j > 0; --j ) {
      for( int i = 0; i < j; ++i ) {
        ++pairs;
        planeStatPerPair = 0;
        switch( algo ) {
          case  0 : if( ! tetDisjointSAT( tets[i], tets[j]) ) ++gNumHits; break;
          case  1 : if( ! tetDisjoint( tets[i], tets[j]) ) ++gNumHits; break;
          case  2 : if( ! sphericalDisjoint( tets[i], tets[j]) ) ++gNumHits; break;
          case  3 : if( ! gjkDisjoint( tets[i], tets[j]) ) ++gNumHits; break;
          case  4 : if( 0.0f == sphericalDistance( tets[i], tets[j]) ) ++gNumHits; break;
          case  5 : if( 0.0f == gjkDistance( tets[i], tets[j]) ) ++gNumHits; break;
          default : if( ! faceSeparation( tets[i], tets[j]) ) ++gNumHits; break;
        }
        planeStat.addSample(planeStatPerPair);
        /*case  4 : if( ! tetDisjointDirect( tets[i], tets1[j]) ) ++gNumHits; break; 
          case  5 : {
          bool overlayHit = ! tetDisjoint( tets[i], tets1[j]);
          bool satHit = ! tetDisjointSAT( tets[i], tets1[j]);
#ifdef WITH_CGAL
bool cgalHit = ! tetDisjointCGAL( tets[i], tets1[j]);
#endif
if( overlayHit != satHit
#ifdef WITH_CGAL
|| overlayHit != cgalHit 
#endif
) {
++gNumHits;
#ifndef WITH_CGAL
cout
<< "Overlay: " << overlayHit
<< ", SAT: " << satHit
        //<< ", CGAL: " << cgalHit
        << endl;
#else
if( overlayHit == cgalHit )
++disRight;
else if( satHit == cgalHit )
++satRight;
else
++bothWrong;
printMistakes();
#endif
} else
++bothRight;
break;
}
*/
        }
}
myChrono.stop();
#ifdef WITH_CGAL
printMistakes();
#endif
cout << gNumHits << " hits over " << pairs << " pairs (";
cout << (100.0*(double)gNumHits)/pairs << "%) in " << myChrono.elapsed_time() << " seconds ( ";
cout << (pairs / myChrono.elapsed_time() / 1000.0) << " KPairs/s ).";
cout << endl << planeStat.min() << " / " << planeStat.mean() << "+/-" << planeStat.dev() << " / " << planeStat.second_max() << " / " << planeStat.max()
    << " min/avg/dev/premax/max directions/pair.";
cout << endl << "FAILS:" << nbFails << endl << flush;
}
closeDatabase();
return EXIT_SUCCESS;
}
