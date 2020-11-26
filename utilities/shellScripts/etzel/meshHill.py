#!/usr/bin/env python

import math
import os

readinVars=['Nx', 'Ny', 'Nz', 'yStar', 'incinterval', 'height', 'dyOnePlus'];
Nmax = 500000;

dict='constant/polyMesh/blockMeshDict'
include='constant/polyMesh/meshParameters'

fDict = open(dict, 'r')

print "Reading dict: ", dict
for line in fDict.readlines():
    data = [x.strip() for x in line.split(None)]
    if (not data):
        continue

    for v in readinVars:
        if data[0] == v:
            if len(data) > 1:
                number = data[1]
                print v, ": ",  number
                # Remove anything after semicolon
                vars()[v] =  number.split(';')[0]
            else:
                print data[0], " not defined!"

print "Finished reading dict.\n"
fDict.close()  # Close dict

# Treat default values
#- Assign dyOnePlus as 1 if not defined in blockMeshDict
if(not ('dyOnePlus' in globals())):
    print "dyOnePlus not defined. Taking yPlus = 1 as default."
    dyOnePlus = 1;

# convert string to numbers
Nx=int(Nx)
Ny=int(Ny)
Nz=int(Nz)
yStar=float(yStar)
#incinterval=float(incinterval)
dyOnePlus=float(dyOnePlus)
incinterval=1.555E-6
height=float(height)
middle=height/2
splitX=(height+1)/2

r = 1.0+1.0e-7;
ni = 1;
while ((ni < Nmax) & (((middle-2.0*dyOnePlus*yStar*(1-(r**(Ny/2)))/(1-r))**2)
       > ((middle-2.0*dyOnePlus*yStar*(1-((r+incinterval)**(Ny/2)))/(1-(r+incinterval)))**2))):
    r=r+0.00001;
    ni += 1;

print "Iterations: ", ni
if(ni == Nmax):
    print "Warning: Iteration quit after ", Nmax, "Iterations"

print "Stretch ratio =", r

Rb = (r**(Ny/2-1));
Rt = 1.0/Rb;

# compute and write other parameters
Nxc= Nx / 3
Nxh= (Nx - Nxc) / 2
Nyt=Ny / 2
Nyb=Ny - Nyt

# compute dx+, dy+, dz+
dxPlus = (9.0/Nx)/yStar;
dzPlus = (4.5/Nz)/yStar;
dyPlusMin = (middle*(1-r)/(1-(r**Nyb)))/yStar
dyPlusMax = dyPlusMin * Rb


fInclude = open(include, 'w')
fInclude.write('Nxc ' + str(Nxc) +';\n');
fInclude.write('Nxh ' + str(Nxh) +';\n');
fInclude.write('Nyt ' + str(Nyt) +';\n');
fInclude.write('Nyb ' + str(Nyb) +';\n');
fInclude.write("middle " + str(middle) + ";\n");
fInclude.write("split " + str(splitX) + ";\n");
fInclude.write("Rt " + str(Rt) + ";\n");
fInclude.write("Rb " + str(Rb) + ";\n");

fInclude.write("\n//Stretch ratio: " + str(r) + "\n")
fInclude.write("//Resolution Information:\n//Dx+ * Dy+ * Dz+ = " \
               + str(dxPlus) + ' x [' + str(dyPlusMin) + ', ' + str(dyPlusMax) \
               +  '] x ' + str(dzPlus) + "\n")
fInclude.write("//Note: dy+ is half of Dy+, that is %.1f \n" % (dyPlusMin*0.5) )
fInclude.write("//Based on smallest wall unit. The number are divided by 2 if based on upper wall unit.")

fInclude.close();

if ( r > 1.05 ):
    print "\nWarning: r>1.05, please increase Ny!"

# print "top section ( 2 < y < 3 ) of channel has ", Nyt, " cells in wall normal direction"
# print "bottom section ( 0 < y < 2 ) of channel has", Nyb, " cells in wall normal direction"
print "Nx * Ny * Nz = ", Nx, " x ",  Ny, " x ",  Nz, "[", Nx*Ny*Nz, "cells]"
print "Dx+ * Dy+ * Dz+ = %.1f x [%.1f, %.1f] x %.1f" % (dxPlus, dyPlusMin, dyPlusMax, dzPlus)
print "Note: dy+ is half of Dy+, that is %.1f" % (dyPlusMin*0.5)
print "Based on smallest wall unit. The number are divided by 2 if based on upper wall unit."
print "Computed mesh parameters written to ", include

os.system('blockMesh > blockMesh.log')
