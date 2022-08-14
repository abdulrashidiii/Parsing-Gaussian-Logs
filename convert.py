#!/usr/bin/env python

"""
Usage: python3 convert.py <log filename> <amino acid name>
"""

from sys import argv, exit
import numpy as np

fname = argv[1]
aa = argv[2].title()

aa_code = {
    'Glycine': ['Gly', 'G', np.array([19,16,14,15,17,18,1,3,4,7,6,2,5,8,10,13,9,11,12])], 
    'Alanine': ['Ala', 'A', np.array([22,19,17,18,20,21,1,3,6,10,4,2,5,7,8,9,11,13,14,12,15,16])], 
    'Valine': ['Val', 'V', np.array([26,25,23,24,27,28,1,3,6,10,4,2,5,8,11,9,14,7,12,13,15,16,17,19,20,18,21,22])], 
    'Cysteine': ['Cys', 'C', np.array([23,20,18,19,21,22,1,3,6,10,4,2,5,9,7,8,11,12,14,15,13,16,17])], 
    'Serine': ['Ser', 'S', np.array([22,19,18,23,20,21,1,3,6,10,4,2,5,9,7,8,11,12,14,15,13,16,17])],
    'Aspartate': ['Asp', 'D', np.array([20,19,18,23,21,22,1,3,6,10,4,2,5,9,11,24,7,8,25,12,14,15,13,16,17])],
    'Aspargine': ['Asn', 'N', np.array([25,22,21,26,23,24,1,3,6,10,4,2,5,9,11,12,13,7,8,14,15,17,20,16,18,19])],
    'Threonine': ['Thr', 'T', np.array([23,22,21,26,24,25,1,3,6,10,4,2,5,8,11,9,7,12,13,14,15,17,18,16,19,20])],
    'Phenylalanine': ['Phe', 'F', np.array([29,28,27,32,30,31,1,3,6,10,4,2,5,9,12,15,18,14,11,7,8,16,19,20,17,13,21,23,24,22,25,26])],
    'Histidine': ['His', 'H', np.array([28,25,24,29,26,27,1,3,6,10,4,2,5,9,12,14,13,11,7,8,15,17,16,18,20,21,19,22,23])],
    'Tyrosine': ['Tyr', 'Y', np.array([30,29,28,33,31,32,1,3,6,10,4,2,5,9,12,15,18,14,11,7,8,16,19,20,21,17,13,22,24,25,23,26,27])],
    'Isoleucine': ['Ile', 'I', np.array([30,27,26,31,29,30,1,3,6,10,4,2,5,9,16,17,8,11,12,13,7,14,15,18,19,20,22,23,21,24,25])],
    'Leucine': ['Leu', 'L', np.array([30,27,26,31,28,29,1,3,6,10,4,2,5,9,12,14,13,17,7,8,11,15,16,18,19,20,22,25,21,23,24])],
    'Tryptophan': ['Trp', 'W', np.array([33,32,31,36,34,35,1,3,6,10,4,2,5,9,11,14,15,18,22,19,16,12,7,8,13,17,21,24,23,20,25,27,28,26,29,30])],
    'Glutamate': ['Glu', 'E', np.array([23,22,21,26,24,25,1,3,6,10,4,2,5,9,13,14,27,7,8,11,12,28,15,17,18,16,19,20])],
    'Methionine': ['Met', 'M', np.array([28,25,24,29,26,27,1,3,6,10,4,2,5,9,13,14,15,7,8,11,12,16,17,18,20,21,19,22,23])],
    'Glutamine': ['Gln', 'Q', np.array([26,25,24,29,27,28,1,3,6,10,4,2,5,9,13,14,15,16,7,8,11,12,17, 18,20,21,19,22,23])],
    'Lysine': ['Lys', 'K', np.array([30,29,28,33,31,32,1,3,6,10,4,2,5,9,13,16,19,20,7,8,11,12,14,15,17,18,21,22,24,25,23,26,27])],
    'Arginine': ['Arg', 'R', np.array([32,31,30,35,33,34,1,3,6,10,4,2,5,9,13,16,18,20,22,19,21,7,8,11,12,14,15,17,23,24,26,27,25,28,29])],
    'Selenocysteine': ['Sec', 'X', np.array([20,19,18,23,21,22,1,3,6,9,4,2,5,11,7,8,10,12,14,15,13,16,17])],
    'Proline': ['Pro', 'P']
    }

try:
  three_letter, one_letter, trans_mat = aa_code[aa][0], aa_code[aa][1], aa_code[aa][2]

except KeyError:
  print('Amino acid name -- %s -- not valid! Please try again and use complete amino acid name.' % aa)
  exit(1)

def distance(x,y):
  return np.linalg.norm(x-y)

def angle(x,y,z):
  xy = x - y 
  zy = z - y

  xy_zy = np.dot(xy, zy)
  xy_xy = np.linalg.norm(xy)
  zy_zy = np.linalg.norm(zy)

  cos_theta = xy_zy / (xy_xy * zy_zy)

  return np.arccos(cos_theta) * 180 / np.pi

def dihedral(x,y,z,a):
  yx = y - x
  zy = y - z
  az = a - z

  n1 = np.cross(yx,zy)
  n1 = n1 / np.linalg.norm(n1)
  n2 = np.cross(zy,az)
  n2 = n2 / np.linalg.norm(n2)

  zy = zy / np.linalg.norm(zy)
  m1 = np.cross(n1,zy)

  x = np.dot(n1,n2)
  y = np.dot(m1,n2)

  return np.arctan2(y,x) * 180 / np.pi

xyz_coords = np.loadtxt(fname, usecols=[1,2,3])
n_atoms = xyz_coords.shape[0]

try:
  assert n_atoms == len(trans_mat)
except AssertionError:
  print('Number of atoms in xyz file do not match expected number. Please recheck.')
  exit(1)

data = np.empty((n_atoms,3))
for i in range(len(xyz_coords)):
  data[i] = xyz_coords[trans_mat[i]-1].reshape(1,3)

job = "# b3lyp/6-31g(d) opt=(modredundant,calcfc,tight)\n"

with open('temp/%s_zmat.txt' % one_letter, 'r') as f:
  zmat = f.readlines()
  try:
    assert n_atoms == len(zmat) - 3
  except AssertionError:
    print('Number of atoms do not match excpected number. Please recheck.')
    exit(1)

  R2 = zmat[4].split()
  R2 = np.array([R2[2][1:], R2[1]], dtype='int') - 1

  A3 = zmat[5].split()
  A3 = np.array([A3[2][1:], A3[1], A3[3]], dtype='int') - 1

  D_mat = np.empty((n_atoms-3, 4), dtype='int')
  D_zmat = zmat[6:]
  for i in np.arange(n_atoms-3):
    D = D_zmat[i].split()
    D_mat[i] = np.array([D[2][1:], D[1], D[3], D[5]], dtype='int') - 1

with open('std_zmat_files/%s_par_opt.txt' % fname[10:-4], 'w') as f:
  f.write('R2 %.5f\n' % distance(data[R2[0]], data[R2[1]]))
  f.write('R3 %.5f\n' % distance(data[A3[0]], data[A3[1]]))
  
  for r in range(4,n_atoms+1):
    i = r - 4
    f.write('R%d %.5f\n' % (r, distance(data[D_mat[i][0]], data[D_mat[i][1]])))
  
  f.write('A3 %.5f\n' % angle(data[A3[0]], data[A3[1]], data[A3[2]]))
  
  for a in range(4,n_atoms+1):
    i = a - 4
    f.write('A%d %.5f\n' % (a, angle(data[D_mat[i][0]], data[D_mat[i][1]], data[D_mat[i][2]])))
  
  for d in range(4,n_atoms+1):
    i = d - 4
    f.write('D%d %.5f\n' % (d, dihedral(data[D_mat[i][0]], data[D_mat[i][1]], data[D_mat[i][2]], data[D_mat[i][3]])))
