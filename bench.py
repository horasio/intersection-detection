#!/usr/bin/env python
# encoding: utf-8

from subprocess import call
import sys, time
from os import path
from random import randint

nbArgs = len(sys.argv)

if nbArgs < 2:
  print("Usage:\n\tbench.py  tet-directory  [tetra,obb,few,frustum,gen50]*")
  exit(-1)

prog = path.join(sys.argv[1], "tet")

print("Benchmarking "+prog)

dbName = time.strftime("statistics-%Y-%m-%d.db", time.localtime())

baseCommand = [prog, "-b", dbName] # -b: write stats to DB. -f: fullscreen

def do(command):
  ## first, simulate a mouse move so that Ubuntu does not freeze the world after several hours of operations.
  #call(["xdotool","mousemove",str(randint(1,100)), str(randint(1,100))])
  #print("-------------------------------------------------------------------------")
  #print(' '.join(command))
  call(command)

nbo = 1000

algorithms = ["naive", "sat", "sphere", "gjk"]
benchmarks = ["tetra", "obb", "few", "frustum", "gen50"]

if nbArgs > 2:
	benchmarks = [sys.argv[i] for i in range(2,nbArgs)]

# TETRAHEDRA

def benchTetrahedra(iternum):
    nbo = 2000
    spreads = [0, 1.68] + [pow(2, i/2.0) for i in range(25)]
    for algo in algorithms:
        for spread in spreads:
            print("TETS [{0}/{1}] ".format(iternum,nbIter) + algo + ' ' + str(spread))
            command = baseCommand + ['--tet', '--number-of-convexes', str(nbo), '--algorithm', algo, '--spread', str(spread)]
            do(command)

# A FEW VERTICES

def benchFew(iternum):
    spreads = [0, 2.4, 3.36] + [pow(2, i/2.0) for i in range(2,26)]
    for algo in ['sphere', 'gjk', 'naive']:
        if algo == 'naive':
            nbo = 300
        else:
            nbo = 1600
        for spread in spreads:
            print("FEW [{0}/{1}] ".format(iternum,nbIter) + algo + ' ' + str(spread))
            command = baseCommand + ['--general', '-v', '16', '--number-of-convexes', str(nbo), '--algorithm', algo, '--spread', str(spread)]
            do(command)

# FRUSTUM CULLING

def benchFrustumCulling(iternum):
    nbo = 1000
    spreads = [0] + [pow(2, i/3.0)/10.0 for i in range(25)]
    spreads = [pow(1.15, i) for i in range(24)]
    for algo in ['greene', 'sphere', 'gjk', 'hybrid']:
        for spread in spreads:
            print("FRUSTUM-AABB [{0}/{1}] ".format(iternum,nbIter) + algo + ' ' + str(spread))
            do(baseCommand + ['--frustum-aabb', '--number-of-convexes', str(nbo), '--algorithm', algo, '--spread', str(spread)])
            if algo == 'greene': continue
            print("FRUSTUM-SPHERE [{0}/{1}] ".format(iternum,nbIter) + algo + ' ' + str(spread))
            do(baseCommand + ['--frustum-sphere', '--number-of-convexes', str(nbo), '--algorithm', algo, '--spread', str(spread)])

# ORIENTED BOUNDING BOXES

def benchOBB(iternum):
    spreads = [0] + [0.5 * pow(2, i/4.0) for i in range(26)]
    nbo = 2000
    for algo in algorithms:
        if algo == "naive":
            continue
        for spread in spreads:
            print("OBB [{0}/{1}] ".format(iternum,nbIter) + algo + ' ' + str(spread))
            command = baseCommand + ['--obb', '--number-of-convexes', str(nbo), '--algorithm', algo, '--spread', str(spread)]
            do(command)

# GENERAL at 50%

def benchGeneralAt50(iternum):
    sizes = [5, 10, 20, 50, 100, 200, 500, 1000, 2000]
    spread = {4:1.4, 5:2.62, 8:4.3, 10:4.8, 16:5.55, 20:5.81, 50:6.4, 100:6.6, 200:6.7, 500:6.8, 1000:6.8, 2000:6.8}
    for size in sizes:
        for algo in ['gjk','sphere','naive']:
            localNbo = 2000
            if algo == "naive":
                if size > 1000: continue
                if size <= 10: localNbo = 2000
                elif size <= 20: localNbo = 800
                elif size <= 50: localNbo = 300
                elif size <= 100: localNbo = 150
                elif size <= 200: localNbo = 80
                elif size <= 500: localNbo = 50
                elif size <= 1000: localNbo = 30
                else: localNbo = 30
            else:
                if size >= 1000: localNbo = 1000
            print("GENERAL [{0}/{1}] ".format(iternum,nbIter) + algo + ' spread:' + str(spread[size]) + ', size:' + str(size))
            command = baseCommand + ['--general', '--number-of-convexes', str(localNbo), '--number-of-vertices', str(size), '--algorithm', algo]
            command = command + ['--spread', str(spread[size])]
            do(command)
            if algo != "naive":
                do(command + ['--sorted-vertices'])

metaIters = 40
iters = 5
nbIter = metaIters * iters
for m in range(metaIters):
	base = 1 + m * iters
	if "tetra" in benchmarks:
		for i in range(iters):
			benchTetrahedra(base+i)
	if "obb" in benchmarks:
		for i in range(iters):
			benchOBB(base+i)
	if "few" in benchmarks:
		for i in range(iters):
			benchFew(base+i)
	if "frustum" in benchmarks:
		for i in range(iters):
			benchFrustumCulling(base+i)
	if "gen50" in benchmarks:
		for i in range(iters):
			benchGeneralAt50(base+i)

