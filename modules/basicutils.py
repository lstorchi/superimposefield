import subprocess
import numpy 
import pybel
import sets
import math
import sys
import re
import os

CONVERTER = 1.889716
DELTAVAL = 10.0
STEPVAL = 1.0

###############################################################################

def bufcount(filename):
    f = open(filename)                  
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization

    buf = read_f(buf_size)
    while buf:
      lines += buf.count('\n')
      buf = read_f(buf_size)

    f.close()

    return lines

###############################################################################

def readkontfile (kontname, probe):

  lineamnt = bufcount(kontname)
  
  dim = (lineamnt - 1)/2

  energy = numpy.empty([1,1,1], float)

  if ((dim * 2) + 1) != lineamnt :
    print "Maybe invalid kont file"
    exit(1)

  fk = open(kontname)

  xsets = sets.Set()
  ysets = sets.Set()
  zsets = sets.Set()
  switchtofieldm = False

  nx = ny = nz = 0
  ix = iy = iz = 0
  for l in fk:

    if (l.find(probe) >= 0):
      switchtofieldm = True 
      nx = len(xsets)
      ny = len(ysets)
      nz = len(zsets)
      energy = numpy.arange(nx*ny*nz, dtype=float).reshape(nx, ny, nz)

    else:
      if switchtofieldm:
        p = re.compile(r'\s+')
        line = p.sub(' ', l)
        line = line.lstrip()
        line = line.rstrip()

        e = float(line)
        energy[ix, iy, iz] = e

        iy = iy + 1
        if (iy == ny):
          iy = 0
          ix = ix + 1
        
        if (ix == nx):
          ix = 0
          iy = 0
          iz = iz + 1

        if (iz == nz):
          ix = 0
          iy = 0
          iz = 0

      else:
        p = re.compile(r'\s+')
        line = p.sub(' ', l)
        line = line.lstrip()
        line = line.rstrip()
        n, x, y, z = line.split(" ")

        xsets.add(float(x))
        ysets.add(float(y))
        zsets.add(float(z))

  fk.close()

  return energy

###############################################################################

def write_to_cube (mol1, mol1field, fname, xnstep, ynstep, znstep,\
    step, xmin, ymin, zmin):

  opf = open(fname, "w")

  print "Writing final cube... "

  zero = 0.0
  opf.write("El field\n")
  opf.write("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (len(mol1.atoms), CONVERTER*xmin, \
      CONVERTER*ymin, CONVERTER*zmin))
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (xnstep, CONVERTER*step, zero, zero))
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (ynstep, zero, CONVERTER*step, zero))
  opf.write("%4d %11.6f %11.6f %11.6f\n" % (znstep, zero, zero, CONVERTER*step))

  for atom in mol1:
    x = CONVERTER*atom.coords[0]
    y = CONVERTER*atom.coords[1]
    z = CONVERTER*atom.coords[2]
    c = atom.partialcharge
    anu = atom.atomicnum

    opf.write("%4d %11.6f %11.6f %11.6f %11.6f\n" % (anu, c, x, y, z))

  for ix in range(xnstep):
    for iy in range(ynstep):
      for iz in range(znstep):
        opf.write("%g "%mol1field[ix,iy,iz])
        if (iz % 6 == 5):
          opf.write("\n")
      opf.write("\n")

  opf.close()


###############################################################################

def get_mean_gridfield (filename, weightfile, probe, rootdir):

   xmin = float("inf")
   ymin = float("inf")
   zmin = float("inf")
   xmax = float("-inf")
   ymax = float("-inf")
   zmax = float("-inf")
   
   mol1list = list(pybel.readfile("mol2", filename))
   
   for conf1 in mol1list:
     for a in conf1.atoms:
       x, y, z = a.coords
   
       if x < xmin:
         xmin = x
       if y < ymin:
         ymin = y
       if z < zmin:
         zmin = z
   
       if x > xmax:
         xmax = x
       if y > ymax:
         ymax = y
       if z > zmax:
         zmax = z
   
   xmin = xmin - DELTAVAL
   xmax = xmax + DELTAVAL
   
   ymin = ymin - DELTAVAL
   ymax = ymax + DELTAVAL
   
   zmin = zmin - DELTAVAL
   zmax = zmax + DELTAVAL
   
   print "Grid used: ", xmin, ymin, zmin, xmax, ymax, zmax
   
   xnsteps = int((xmax - xmin) / 1.0) + 1
   
   weights = []
   weightsfp = open(weightfile)
   weightsfp.readline()
   idx = 0
   sum = 0.0
   for line in weightsfp:
     values = re.split('\s+', line)
     weights.append(float(values[6]))
   
     sum = sum + float(values[6])
   
     idx = idx + 1
   
   for i in range(len(weights)):
     weights[i] = weights[i] / sum
   
   weightsfp.close()
   
   if (len(mol1list) != len(weights)):
     print "Dimension error ", len(mol1list) , " vs " , \
       len(weights)
     exit(1)
   
   if rootdir != "./":
       subprocess.call("ln -s "+rootdir+"fixpdb ./", shell=True)
       subprocess.call("ln -s "+rootdir+"kouttokont ./", shell=True)

   energy = numpy.empty([1,1,1], float)
   globalindex = 0
   for conf1 in mol1list:
     outfl = filename
     outfl = outfl.replace(".mol2", "_")
     basicfname = outfl+str(globalindex)+"_"+probe
   
     output = pybel.Outputfile("pdb", basicfname+".pdb")
     output.write(conf1)
     output.close()
   
     toexe = "./fixpdb --remove-all-H2O --unkn-residue-to-grid-types --kout-out="+ \
         basicfname+".kout "+basicfname+".pdb"
     print toexe
     subprocess.call(toexe, shell=True)
   
     kontname = basicfname+".kont"
   
     limits = str(xmin)+";"+str(ymin)+";"+str(zmin)+";"+ \
             str(xmax)+";"+str(ymax)+";"+str(zmax)
     torun = "./kouttokont -g \""+limits+"\" "+basicfname+".kout "+probe + \
             " > "+kontname
     print torun
     subprocess.call(torun, shell=True)
   
     os.remove(basicfname+".pdb")
     os.remove(basicfname+".kout")
   
     # read kont file
     energy1 = readkontfile(kontname, probe)
   
     print "nx: ", energy1.shape[0], " ny: ", energy1.shape[1], \
         " nz: ", energy1.shape[2]
   
     os.remove(kontname)
   
     print "Dealing with: ", kontname, " w: ", weights[globalindex]
   
     if  globalindex == 0:
       nx = energy1.shape[0]
       ny = energy1.shape[1]
       nz = energy1.shape[2]
       energy = numpy.arange(nx*ny*nz, dtype=float).reshape(nx, ny, nz)
       energy = numpy.zeros([nx,ny,nz], float)
   
     energy = energy + weights[globalindex] * energy1
   
     globalindex = globalindex + 1

   return energy, [xmin, ymin, zmin, xmax, ymax, zmax]
