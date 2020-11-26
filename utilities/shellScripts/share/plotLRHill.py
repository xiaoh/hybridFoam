#!/usr/bin/env python

import sys
import os.path as osp
import os
import numpy as np
import matplotlib
#matplotlib.use('PDF')
import matplotlib.pyplot as plt


if len(sys.argv) < 2 :
    datfile = 'CaseNames.in'
else :
    datfile = sys.argv[1]

projectName, dummy = osp.splitext(datfile)
print "project name: ", projectName
f = open(datfile, 'r')

recs = f.read().split(';')

datDir = recs[0]

figdir = osp.join(os.getcwd(), 'figs-' + projectName)
if(not osp.isdir(figdir)):
    os.mkdir(figdir)

allCases = []
allLabels = []
allLines = []
for i in range(1, (len(recs))):
    # Split to 3-element records: Casename, legend, line type   
    rec=recs[i].split()
    if(len(rec)==3):
        case, legend, line = rec
        allCases.append(case)
        allLabels.append(legend)
        allLines.append(line)


colors = ['black', 'blue', 'red']
widths = [2, 2 , 2, 2]
lines  = ['solid', 'solid', 'solid', 'dashdot']

allLocs = map(str, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
allComps = ['U', 'uu', 'uv', 'k']
xLabels = ['$U/U_b$', '$<uu>/U_b^2$', '$<uv>/U_b^2$', '$k/U_b^2$']
allCols = [1, 3, 5, 6]

fignum = []
figname = []
for i in range(len(allCases)):
    for k in range(len(allLocs)):
        loc = allLocs[k]
        #xyfile = 'x'+str(loc)+'_normalised_Results.xy'
        xyfile = 'x'+str(loc)+'_normalized_Results.xy'
        case = osp.join(datDir, allCases[i], xyfile)
        f=open(case, 'r') 
        # Load data
        dataRaw = np.loadtxt(f, comments='%')
        data = np.array(dataRaw)
    
        for j in range(len(allComps)):
            fn = k * len(allComps) + j
            fignum.append(fn)
            
            col = allCols[j]
            comp = allComps[j]
            figname.append(comp+'-'+loc)
            plt.figure(fn)
            plt.plot(data[:, col], data[:, 0], allLines[i], label=allLabels[i], lw=2)
            plt.title(projectName + ": " + xLabels[j] + ' @ $x =$ ' + loc)
            plt.xlabel(xLabels[j])
            plt.ylabel('$y$')
            ymin = data[:, 0].min()
            ymax = data[:,0].max()
            #plt.ylim([ymin, ymax])


# draw legend, title, and then save figures
for i in (range(len(fignum))):
    plt.figure(fignum[i])
    plt.legend(loc=0)
    plt.savefig(osp.join(figdir, figname[i] +'.png'))

# show figures    
# plt.show()





