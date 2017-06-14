#!/usr/bin/env python3.4
# encoding: utf-8

import sys, time

import sqlite3
import math
import numpy
import matplotlib
#matplotlib.use('module://backend_ipe')
import matplotlib.pyplot as plt

genPNG = True

nbArgs = len(sys.argv)

if nbArgs < 2:
  print("Usage:\n\tmakeplot.py  database  [tetra,obb,few,fa,fs,gen50,fg]*")
  exit(-1)

dbName = sys.argv[1]

shortDate = dbName[-8:-3]
tableName = 'statistics'

benchmarks = ["tetra", "obb", "few", "fa", "fs", "gen50", "fg"]
if nbArgs > 2:
	benchmarks = [sys.argv[i] for i in range(2,nbArgs)]

print("Plotting "+dbName)

db = sqlite3.connect(dbName)
cursor = db.cursor()

dpiValue = 300

# =============================================================================

def constraintsToSQL(constraints):
    cons = []
    for key, val in constraints.items():
        if key[0] == '!':
            cons.append('<>'.join([key[1:], str(val)]))
        else:
            cons.append('='.join([key, str(val)]))
    return ' AND '.join(cons)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def maxLabel(t):
    if t == 'gjk':
        return "max"#=MAX_ITER"
    else:
        return "max"

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def title(test_name, t):
    dico = {
            'sat':'SAT',
            'naive':'Generic',
            'sphere':'DSS',
            'sphere1':'DSS, sorted',
            'gjk':'GJK',
            'gjk1':'GJK, sorted',
            'hybrid':'Hybrid',
            'greene':'Greene'
            }
    if test_name == "'TETS'":
        dico['naive'] = 'Tetra'
    try:
        return dico[t]
    except:
        return t

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def color(t):
    dico = {
            'sat':'b',
            'naive':'g',
            'sphere':'r',
            'sphere1':'#ff5555',
            'gjk':'purple',
            'gjk1':'#ff55ff',
            'greene':'#55dd55',
            'hybrid':'y'
            }
    try:
        return dico[t]
    except:
        return 'black'

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def dev_color(t):
    dico = {'sat'   :(0.0, 0.0, 1.0, 0.2),
            'naive' :(0.0, 1.0, 0.0, 0.2),
            'sphere':(1.0, 0.0, 0.0, 0.2),
            'gjk'   :(1.0, 0.0, 1.0, 0.2),
            'hybrid':(1.0, 1.0, 0.0, 0.2)}
    try:
        return dico[t]
    except:
        return 'black'


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def filtered(l,pos):
    vs = [v[pos] for v in l]
    lenl = len(l)
    if lenl == 0:
         print('Empty')
         return l
    avg = sum(vs)/lenl
    n = []
    thres = 1.5
    invt = 1.0 / thres
    for v in l:
        ratio = v[pos] / avg
        if ratio < thres and ratio > invt:
            n.append(v)
    lenn = len(n)
    if lenn == 0:
        print("Too noisy, NOT filtered !")
        return l
    if lenn < lenl:
        print(lenl, "->", lenn)
    else:
        print(lenl, "->  OK")
    return n

def mean(l):
    n = float(len(l))
    return sum(l)/n

def dev(l,m):
    n = float(len(l))
    sqmean = functools.reduce(lambda a, b: x + y*y, l, 0) / n
    # 'm' is assumed to be the mean of the same input list
    return math.sqrt(sqmean - (m * m))

# =============================================================================

def plotSpeedVSVertices(constraints, algos, caption, filename, ymin=900000.0, ymax=40000000.0, linear=False, xlabel=None):
    test_name = constraints['TEST_NAME']
    commandTmpl = "select distinct N_VERTEX from {0} WHERE {1} ORDER BY N_VERTEX ASC"
    command = commandTmpl.format(tableName, constraintsToSQL(constraints))
    cursor.execute(command)
    nbvs = [ s[0] for s in cursor ]
    commandTmpl = "select N_HITS, N_PAIR, TIME from {0} WHERE {1}".format(tableName, "{0}")
    plot_data = {}
    for algo in algos:
        plot_data[algo] = ([], [])
        data = plot_data[algo]
        if algo[-1] == '1':
            algo = algo[:-1]
            constraints['SORTED'] = 1
        else:
            constraints['SORTED'] = 0
        constraints['ALGORITHM'] = "'"+algo+"'"
        for nbv in nbvs:
            constraints['N_VERTEX'] = nbv
            command = commandTmpl.format(constraintsToSQL(constraints))
            cursor.execute(command)
            tot_hits = 0
            tot_pairs = 0
            tot_time = 0.0
            values = filtered(list(cursor), 2)
            for c in values:
                tot_hits += c[0]
                tot_pairs += c[1]
                tot_time += c[2]
            if tot_time > 0.0:
                data[0].append(nbv)#tot_hits/float(tot_pairs))
                data[1].append(tot_pairs / tot_time)
    if 'sphere1' in algos:
        sph1 = plot_data['sphere1'][1]
        gjk1 = plot_data['gjk1'][1]
        sph = plot_data['sphere'][1]
        gjk = plot_data['gjk'][1]
        print("max sphere/gjk ratio=", max([sph[s]/gjk[s] for s in range(len(sph))]))
        print([(s,sph[s]/gjk[s]) for s in range(len(sph))])
        print([(s,sph1[s]/gjk1[s]) for s in range(len(sph))])
        print([(s,v) for s,v in enumerate(plot_data['sphere'][0])])
    fig, ax = plt.subplots(figsize=(10,4), dpi=dpiValue)
    ax.set_xscale('log')
    if linear:
        ax.set_yscale('linear')
    else:
        ax.set_yscale('log')
    ax.grid(b=True, which='major', color='black')#, linestyle='-')
    ax.grid(b=True, which='minor', color='gray')#, linestyle='-')
    if xlabel == None:
	    ax.set_xlabel('Number of vertices')
    else:
	    ax.set_xlabel(xlabel)
    ax.set_ylabel('Tested pairs per second')
    ax.set_ylim([ymin, ymax])
    ax.set_xlim([4, 2500])
    if xlabel == None:
	    plt.xticks(nbvs, [str(n) for n in nbvs])
    else:
	    plt.xticks(nbvs,[])
    lines = []
    names = []
    for algo in algos:
        d = plot_data[algo]
        x_data = d[0]
        #if algo == 'gjk':
        #    x_data = plot_data['sphere'][0]
        lines.append(ax.plot(x_data, d[1], linestyle='-', color=color(algo), linewidth=1.5, marker='.')[0]) # marker='.'
        names.append(title(test_name, algo))
    #x_d = plot_data['naive'][0]
    #y_d = plot_data['sphere'][1]
    #lines.append(ax.plot(x_d, [y_d[0]*x_d[0]/x for x in x_d], linestyle='-', color='yellow', linewidth=1.5, marker='.')[0])
    #
    ax.set_title(caption)
    ax.legend( lines, names, loc='lower left' )
    fig.tight_layout(pad=0.3) # unit is font-size
    if genPNG:
	    plt.savefig(filename.replace('pdf', 'png'), dpi=dpiValue)
    else:
	    plt.savefig(filename, dpi=dpiValue)
    plt.close(fig)

# =============================================================================

def plotSpeedVSDensity(constraints, algos, caption, filename, ymin=900000.0, ymax=40000000.0, linear=False):
    test_name = constraints['TEST_NAME']
    commandTmpl = "select distinct SPREAD from {0} WHERE {1} ORDER BY SPREAD DESC"
    command = commandTmpl.format(tableName, constraintsToSQL(constraints))
    cursor.execute(command)
    spreads = [ s[0] for s in cursor ]
    commandTmpl = "select N_HITS, N_PAIR, TIME from {0} WHERE {1}".format(tableName, "{0}")
    plot_data = {}
    for algo in algos:
        plot_data[algo] = ([], [])
        data = plot_data[algo]
        constraints['ALGORITHM'] = "'"+algo+"'"
        for sp in spreads:
            constraints['SPREAD'] = sp
            command = commandTmpl.format(constraintsToSQL(constraints))
            cursor.execute(command)
            tot_hits = 0
            tot_pairs = 0
            tot_time = 0.0
            values = filtered(list(cursor), 2)
            for c in values:
                tot_hits += c[0]
                tot_pairs += c[1]
                tot_time += c[2]
            if tot_pairs == 0:
                print("No data for algorithm", algo)
                tot_pairs = 1
                tot_time = 1.0
            data[0].append(tot_hits/float(tot_pairs))
            data[1].append(tot_pairs / tot_time)
    if 'sphere' in algos and 'greene' in algos:
        sph_speed = plot_data['sphere'][1]
        greene_speed = plot_data['greene'][1]
        ratios = [greene_speed[s]/sph_speed[s] for s in range(len(sph_speed))]
        print("max Greene/sphere ratio=", max(ratios))
        print("min Greene/sphere ratio=", min(ratios))
        print(ratios)
    if 'sphere' in algos and 'sat' in algos:
        sph_speed = plot_data['sphere'][1]
        sat_speed = plot_data['sat'][1]
        gjk_speed = plot_data['gjk'][1]
        ratiosDSS = [sat_speed[s]/sph_speed[s] for s in range(len(sph_speed))]
        ratiosGJK = [sat_speed[s]/gjk_speed[s] for s in range(len(sph_speed))]
    if 'sphere' in algos and 'gjk' in algos:
        sph_speed = plot_data['sphere'][1]
        gjk_speed = plot_data['gjk'][1]
        ratiosDG  = [sph_speed[s]/gjk_speed[s] for s in range(len(sph_speed))]
        print("max DSS/GJK ratio=", max(ratiosDG))
    fig, ax = plt.subplots(figsize=(10,4), dpi=dpiValue)
    ax.set_xscale('log')
    if linear:
        ax.set_yscale('linear')
    else:
        ax.set_yscale('log')
    ax.grid(b=True, which='major', color='black')#, linestyle='-')
    ax.grid(b=True, which='minor', color='gray')#, linestyle='-')
    ax.set_xlabel('Collision density')
    ax.set_ylabel('Tested pairs per second')
    ax.set_ylim([ymin, ymax])
    ax.set_xlim([0.0005, 1.2])
    plt.xticks([0.001, 0.01, 0.1, 1], [r'0.1 %', r'1 %', r'10 %', r'100 %'])
    lines = []
    names = []
    for algo in algos:
        d = plot_data[algo]
        x_data = d[0]
        #if algo == 'gjk':
        #    x_data = plot_data['sphere'][0]
        lines.append(ax.plot(x_data, d[1], linestyle='-', color=color(algo), linewidth=1.5, marker='.')[0]) # marker='.'
        names.append(title(test_name, algo))
    #
    ax.set_title(caption)
    ax.legend(lines, names, loc='lower left')
    fig.tight_layout(pad=0.3) # unit is font-size
    if genPNG:
	    plt.savefig(filename.replace('pdf', 'png'), dpi=dpiValue)
    else:
	    plt.savefig(filename, dpi=dpiValue)
    plt.close(fig)

# =============================================================================

max_iter = 20

def plotPlaneTests(constraints, algos, caption, filename, ymin, ymax, fdmin, fdmax, nomax=False, nofail=True):
    test_name = constraints['TEST_NAME']
    commandTmpl = "select distinct SPREAD from {0} WHERE {1} ORDER BY SPREAD DESC"
    command = commandTmpl.format(tableName, constraintsToSQL(constraints))
    cursor.execute(command)
    spreads = [ s[0] for s in cursor ]
    commandTmpl = "select N_PLANE, N_SQ_PLANE, N_HITS, N_PAIR, N_PLANE_MAX, N_FAIL from {0} WHERE {1}".format(tableName, "{0}")
    plot_data = {}
    for algo in algos:
        plot_data[algo] = ([], [], [], [], [])
        data = plot_data[algo]
        constraints['ALGORITHM'] = "'"+algo+"'"
        s = []
        for sp in spreads:
            constraints['SPREAD'] = sp
            command = commandTmpl.format(constraintsToSQL(constraints))
            cursor.execute(command)
            tot_planes = 0
            tot_sq_planes = 0
            tot_hits = 0
            tot_pairs = 0
            max_planes = 0
            tot_fail = 0
            values = list(cursor)
            print(algo,":",len(values),"samples. ",end="")
            for c in values:
                tot_planes += c[0]
                tot_sq_planes += c[1]
                tot_hits += c[2]
                tot_pairs += c[3]
                local_max = c[4]
                if (algo not in ['sphere', 'gjk']) or local_max < max_iter:
                    max_planes = max(max_planes, local_max)
                tot_fail += c[5]
            #if max_planes == 0:
            #    max_planes = max_iter-1
            density = tot_hits/float(tot_pairs)
            mean_fail_per_pair = tot_fail/float(tot_pairs)
            mean_plane_tests = tot_planes/float(tot_pairs)
            dev_plane_tests = math.sqrt(tot_pairs * tot_sq_planes - tot_planes * tot_planes) / tot_pairs
            data[0].append(density)
            data[1].append(mean_plane_tests)
            data[2].append(dev_plane_tests)
            data[3].append(max_planes)
            data[4].append(mean_fail_per_pair)
            if tot_fail > 0:
                print(algo, "failures: "+str(tot_fail), "at density", density)
            else:
                print("")
        plot_data[algo] = (data[0], numpy.array(data[1]), numpy.array(data[2]), data[3], data[4])
    fig, ax = plt.subplots(figsize=(10,4), dpi=dpiValue)
    ax.set_xscale('log')
    low_y = ymin
    if False:
        ax.set_yscale('log')
        low_y = 1.0
    ax.grid(b=True, which='major', color='black')#, linestyle='-')
    ax.grid(b=True, which='minor', color='gray')#, linestyle='-')
    ax.set_xlabel('Collision density')
    ax.set_ylabel('Number of plane tests per pair')
    ax.set_ylim([low_y, ymax])
    ax.set_xlim([0.0005, 1.2])
    plt.xticks([0.001, 0.01, 0.1, 1], [r'0.1 %', r'1 %', r'10 %', r'100 %'])
    if not nofail:
        axright = ax.twinx()
        axright.set_xlim([0.0005, 1.2])
        axright.set_ylim([fdmin, fdmax])
        sf = plt.ScalarFormatter(useMathText=True)
        sf.set_scientific(True)
        sf.set_powerlimits((-3,4))
        axright.yaxis.set_major_formatter(sf)
        ytiks = [0.0, 0.00005, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005]
        ytiks = [y for y in ytiks if y <= fdmax]
        axright.set_yticks(ytiks)
        #axright.set_yticklabels(['$0$', '$5\\cdot 10^{-5}$', '0.0001', '$2\\cdot10^{-4}$', '$3\\cdot10^{-4}$', '$5\\cdot10^{-4}$'])
        #axright.set_yticklabels([sf.format_data(v) for v in ytiks])
        axright.set_ylabel('Failure rate ($\\times 10^{-4}$)')
    # plot
    lines = []
    for algo in algos:
         pd = plot_data[algo]
         ax.fill_between(pd[0], numpy.maximum(low_y, pd[1]-pd[2]), pd[1]+pd[2], linestyle='-', color=dev_color(algo), linewidth=0)
         lines.append(ax.plot(pd[0], pd[1], linestyle='-', color=color(algo), linewidth=1.5))
         if not nomax:
              lines.append(ax.plot(pd[0], pd[3], linestyle='--', color=color(algo), linewidth=1.5))
         if not nofail:
              lines.append(axright.plot(pd[0], pd[4], linestyle='-', color=color(algo), linewidth=1.5, marker='o'))
    ax.set_title(caption)
    if nomax:
         legends = sum([[title(test_name, algo), 'fail rate'] for algo in algos], [])
    else:
         legends = sum([[title(test_name, algo), maxLabel(algo)] for algo in algos], [])
    ax.legend( [l[0] for l in lines], legends, loc='center left', ncol=len(algos))
    fig.tight_layout(pad=0.3) # unit is font-size
    if genPNG:
	    plt.savefig(filename.replace('pdf', 'png'), dpi=dpiValue)
    else:
	    plt.savefig(filename, dpi=dpiValue)
    #plt.savefig(filename.replace('pdf', 'ipe'), format="ipe", dpi=dpiValue)
    plt.close(fig)

# =============================================================================

if "tetra" in benchmarks: # TETRAHEDRA
    plotSpeedVSDensity({'TEST_NAME' : "'TETS'"}, ['sat', 'naive', 'gjk', 'sphere'],
                        "Testing pairs of TETRAHEDRA, LOG-LOG diagram",
                        "tet-speed-vs-density.pdf",
                        900000.0, 100000000.0, linear=False)
    plotPlaneTests({'TEST_NAME' : "'TETS'"}, ['sat', 'naive', 'gjk', 'sphere'],
                        "Testing pairs of TETRAHEDRA, LOG-LIN diagram",
                        "tet-plane-tests-vs-density.pdf",
                        0.0, 45.0, 0.0, 0.00015)
    plotPlaneTests({'TEST_NAME' : "'TETS'"}, ['gjk', 'sphere'],
                        "Testing pairs of TETRAHEDRA, LOG-LIN diagram",
                        "tet-gjk-plane-tests-vs-density.pdf",
                        0.0, 4.5, 0.0, 0.00015, nomax=True, nofail=False)
if "obb" in benchmarks: # OBB
    plotSpeedVSDensity({'TEST_NAME' : "'OBB'"}, ["sphere", "gjk", "sat"],
                        "Testing pairs of oriented boxes, LOG-LOG diagram",
                        "obb-speed-vs-density.pdf",
                        2900000.0, 31000000.0, linear=False)
    plotPlaneTests({'TEST_NAME' : "'OBB'"}, ['gjk', 'sphere'],
                        "Testing pairs of oriented boxes, LOG-LIN diagram",
                        "obb-gjk-plane-tests-vs-density.pdf",
                        0.0, 5.0, 0.0, 0.00015, nomax=True, nofail=False)
if "gen50" in benchmarks: # GENERAL at 50% DENSITY
    plotSpeedVSVertices({'TEST_NAME' : "'GENERAL'", '!N_VERTEX': 16}, ["naive", "sphere", "gjk"],#"sphere1", "gjk", "gjk1"],
            "Random polytopes at 50% collision density. LOG-LOG diagram",
            "fifty-speed-vs-vertices.pdf",
            9.0, 7000000.0, linear=False)
    plotSpeedVSVertices({'TEST_NAME' : "'GENERAL'", '!N_VERTEX': 16}, ["sphere", "sphere1", "gjk", "gjk1"],
            "",#"Random polytopes at 50% collision density. LOG-LOG",
            "fifty-speed-vs-vertices-sorted.pdf",
            90000.0, 6000000.0, linear=False, xlabel='')
if "few" in benchmarks: # FEW VERTICES
    plotSpeedVSDensity({'TEST_NAME' : "'GENERAL'", 'N_VERTEX':16}, ["naive", "sphere", "gjk"],
            "Testing pairs of 16-vertices polytopes. LOG-LOG diagram",
            "gen-speed-vs-density.pdf",
            900000.0, 30000000.0, linear=False)
    plotPlaneTests({'TEST_NAME' : "'GENERAL'", 'N_VERTEX':16}, ['sphere', 'gjk'],
            "Testing pairs of 16-vertices polytopes. LOG-LIN diagram",
            "gen-plane-tests-vs-density.pdf",
            0.0, 5.0, 0.0, 0.00015, nomax=True, nofail=False)
if "fs" in benchmarks: # FRUSTUM / SPHERE
    plotSpeedVSDensity({'TEST_NAME' : "'FRUSTUM-SPHERE'"}, ["sphere", "hybrid", "gjk"],
            "Testing pairs of frustum/sphere. LOG-LOG diagram",
            "fs-speed-vs-density.pdf",
            4900000.0, 32000000.0, linear=False)
    plotPlaneTests({'TEST_NAME' : "'FRUSTUM-SPHERE'"}, ['gjk', 'sphere'],
                        "Testing pairs of frustum/sphere, LOG-LIN diagram",
                        "fs-plane-tests-vs-density.pdf",
                        0.0, 4, 0.0, 0.0005, nomax=True, nofail=False)
    plotPlaneTests({'TEST_NAME' : "'FRUSTUM-SPHERE'"}, ['gjk', 'hybrid'],
                        "Testing pairs of frustum/sphere, LOG-LIN diagram",
                        "fs-hybrid-plane-tests-vs-density.pdf",
                        0.0, 4, 0.0, 0.0005, nomax=True, nofail=False)
if "fa" in benchmarks: # FRUSTUM / AAB
    print("============================\nfa\n")
    plotSpeedVSDensity({'TEST_NAME' : "'FRUSTUM-AABB'"}, ["sphere", "hybrid", "gjk"],
            "Testing pairs of frustum/AAB. LOG-LOG diagram",
            "fa-speed-vs-density.pdf",
            5000000.0, 52000000.0, linear=False)
    plotPlaneTests({'TEST_NAME' : "'FRUSTUM-AABB'"}, ['gjk', 'sphere'],
                        "Testing pairs of frustum/AAB, LOG-LIN diagram",
                        "fa-plane-tests-vs-density.pdf",
                        0.0, 4, 0.0, 0.0005, nomax=True, nofail=False)
    plotPlaneTests({'TEST_NAME' : "'FRUSTUM-AABB'"}, ['gjk', 'hybrid'],
                        "Testing pairs of frustum/AAB, LOG-LIN diagram",
                        "fa-hybrid-plane-tests-vs-density.pdf",
                        0.0, 4, 0.0, 0.0005, nomax=True, nofail=False)
if "fg" in benchmarks: # FRUSTUM / AAB
    print("============================\nfg\n")
    plotSpeedVSDensity({'TEST_NAME' : "'FRUSTUM-AABB'"}, ["sphere", "hybrid", "gjk", "greene"],
            "Testing pairs of frustum/AAB. LOG-LOG diagram",
            "fg-speed-vs-density.pdf",
            5000000.0, 210000000.0, linear=False)
