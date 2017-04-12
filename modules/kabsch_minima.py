import numpy.linalg
import numpy
import sys

##############################################################################

def dumpo_matrix (m):

  for i in range(m.shape[0]):
    for j in range(m.shape[1]):
      sys.stdout.write("%10.5f "%(m[i, j]))
    sys.stdout.write("\n")

##############################################################################

def return_rotation_matrix (coord_from=None, coord_to=None, verbose=False):
  a = numpy.zeros((3, 3), numpy.float64)

  centroid_from = coord_from.sum(0) / coord_from.shape[0]
  centroid_to = coord_to.sum(0) / coord_to.shape[0]

  # Loop over the atoms.
  for i in range(coord_from.shape[0]):
    orig_from = coord_from[i] - centroid_from
    orig_to = coord_to[i] - centroid_to
    a += numpy.outer(orig_from, orig_to)
                                          
  u, s, v = numpy.linalg.svd(a)
                                             
  d = numpy.diag([1, 1, numpy.sign(numpy.linalg.det(a))])
                                                        
  r = numpy.dot(numpy.transpose(v), \
          numpy.dot(d, numpy.transpose(u)))

  if (verbose):
    print "Verify rotation matrix "
    dumpo_matrix (r)
    print "det: (should be +/-1)" , numpy.linalg.det(r)
    print "transpose matrix should be equal to inverse matrix" 
    tr = numpy.transpose(r)
    dumpo_matrix (tr)
    print " "
    ir = numpy.linalg.inv(r)
    dumpo_matrix (ir)
            
  return r, centroid_to, centroid_from

##############################################################################

def rmsd(vin, win):

  d = len(vin[0])
  n = len(vin)

  rmsd = 0.0
  for v, w in zip(vin, win):
    rmsd += sum([(v[i]-w[i])**2.0 for i in range(d)])
   
  return numpy.sqrt(rmsd/n)

##############################################################################

def centroid(x):
  
  c = sum(x)/len(x)

  return c
