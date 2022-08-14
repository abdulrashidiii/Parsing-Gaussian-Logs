#!/usr/bin/env python

"""
Usage: python3 blocks.py <log filename>
"""

from sys import argv
import numpy as np

fname = argv[1]

def Distance(x,y):
  return np.linalg.norm(x-y)

def Angle(x,y,z):
  xy = x - y
  zy = z - y

  xy_zy = np.dot(xy,zy)
  xy_xy = np.linalg.norm(xy)
  zy_zy = np.linalg.norm(zy)

  cos_theta = xy_zy / (xy_xy * zy_zy)

  return np.arccos(cos_theta) * 180.0 / np.pi

def Dihedral(x,y,z,a):
  yx = y - x
  zy = y - z
  az = a - z

  n1 = np.cross(yx,zy)
  n1 = n1 / np.linalg.norm(n1)
  n2 = np.cross(zy,az)
  n2 = n2 / np.linalg.norm(n2)

  zy = zy / np.linalg.norm(zy)
  m1 = np.cross(n1, zy)

  x = np.dot(n1, n2)
  y = np.dot(m1, n2)

  return np.arctan2(y,x) * 180.0 / np.pi

def ParseOneOpt(num_opt, lines):
  i = 0
  line = lines[i].strip()
  standard_xyzs = []
  while 'Optimized Parameters' not in line:
    if 'Standard orientation:' not in line:
      i = i + 1
      line = lines[i].strip()
    else:
      standard_xyz = []
      xyz_tmp_lines = lines[i+5:i+5+n_atoms]

      for xyz_tmp_line in xyz_tmp_lines:
        standard_xyz.append(xyz_tmp_line.split()[3:6])

      standard_xyz = np.array(standard_xyz, dtype='float64')
      standard_xyzs.append(standard_xyz)

      i = i + 5 + n_atoms
      line = lines[i].strip()

  tmp_lines = lines[:i]
  tmp_lines2 = lines[i:i+5]
  lines = lines[i+5:]

  tmp = []
  for i in tmp_lines:
    line = i.split()
    if 'SCF' in line:
      tmp.append(line[4])

  energies.append(float(tmp[-1]))

  params_R = []
  params_A = []
  params_D = []
  params_value = []

  i = 0
  line = lines[i].strip()
  while line[0] == '!':
    line = line.strip('!').split()
    temp = line[1][1:].strip('()').split(',')
    if line[0][0] == 'R':
      params_R.append(np.array(temp, dtype='int'))
    elif line[0][0] == 'A':
      params_A.append(np.array(temp, dtype='int'))
    elif line[0][0] == 'D':
      params_D.append(np.array(temp, dtype='int'))
    params_value.append(float(line[2]))
    i = i + 1
    line = lines[i].strip()

  params_R = np.array(params_R) - 1
  params_A = np.array(params_A) - 1
  params_D = np.array(params_D) - 1
  params_value = np.array(params_value, dtype='float64')

  scan_params = []
  for j in range(len(params_D)):
    for k in range(len(modredundant)):
      if np.all(np.equal(params_D[j], modredundant[k])):
        scan_params.append(float(params_value[len(params_R) + len(params_A) + j]))

  g = open('par_files/%d_%s_%s_par_%.2f_%.2f.txt' % (num_opt, title, theory, scan_params[0], scan_params[1]), 'w')
  g.writelines(tmp_lines2)
  g.writelines(lines[:i])
  g.write(lines[i])
  lines = lines[i:]
  g.write("E = %.10f\n" % energies[-1])
  g.close()

  i = 0
  while 'GradGradGrad' not in lines[i]:
    i = i + 1

  lines = lines[i+7:]

  standard_xyzs = standard_xyzs[-1]
  with open('xyz_files/%d_%s_%s_xyz_%.2f_%.2f.xyz' % (num_opt, title, theory, scan_params[0], scan_params[1]), 'w') as xyz_file:
    for atom_name, standard_xyz in list(zip(atoms, standard_xyzs)):
      xyz_file.write('%s\t\t\t%15.10f\t%15.10f\t%15.10f\n' % (atom_name, standard_xyz[0], standard_xyz[1], standard_xyz[2]))
  
  return lines
  

energies = []
with open(fname, 'r') as f:
  lines = f.readlines()

  i = 0
  while '#' not in lines[i].strip():
    i = i + 1

  lines = lines[i:]

  mode = []
  i = 0
  while '------' not in lines[i].strip():
    mode.append(lines[i].strip())
    i = i + 1

  mode = ''.join(mode)
  theory = ''.join(mode.split()[0][1:].split('/'))
  
  lines = lines[i+1:]

  i = 0
  while '-' != lines[i].strip()[0]:
    i = i + 1

  lines = lines[i+1:]

  if 'Title Card Required' not in lines[0].strip():
    title = lines[0].strip()
  else:
    title = fname.split('/')[-1][:-4]

  lines = lines[3:]

  tmp = lines[0].strip().split()

  charge = int(tmp[2])
  multiplicity = int(tmp[5])

  lines = lines[1:]

  i = 0
  line = lines[i].strip()
  atoms = []
  bond_pairs = []
  angle_pairs = []
  dihedral_pairs = []

  if len(lines[0].split()) == 4:
    coords = []
    i = 0
    line = lines[i].split()
    while len(line) != 0:
      atoms.append(line[0])
      coords.append(line[1:4])
      i = i + 1
      line = lines[i].split()
    lines = lines[i+1:]
    coords = np.array(coords, dtype='float64')
    n_atoms = coords.shape[0]
  elif len(lines[0].split()) == 1:
    while line[:9] != 'Variables':
      line = line.split()
      atoms.append(line[0])
      if len(line) > 1:
        bond_pairs.append(int(line[1]))
      if len(line) > 3:
        angle_pairs.append(int(line[3]))
      if len(line) > 5:
        dihedral_pairs.append(int(line[5]))
      i = i+1
      line = lines[i].strip()
  
    n_atoms = i
  
    bond_pairs = np.array(bond_pairs, dtype='int')
    angle_pairs = np.array(angle_pairs, dtype='int')
    dihedral_pairs = np.array(dihedral_pairs, dtype='int')
  
    bond_pairs = np.c_[bond_pairs, (np.arange(n_atoms)+1)[1:]]
  
    angle_pairs = np.c_[angle_pairs, bond_pairs[1:]]
  
    dihedral_pairs = np.c_[dihedral_pairs, angle_pairs[1:]]
  
    R = np.empty(n_atoms - 1, dtype='float64')
    A = np.empty(n_atoms - 2, dtype='float64')
    D = np.empty(n_atoms - 3, dtype='float64')
  
    lines = lines[n_atoms+1:]
  
    for i in range(n_atoms-1):
      R[i] = float(lines[i].strip().split()[1])
  
    lines = lines[n_atoms-1:]
  
    for i in range(n_atoms-2):
      A[i] = float(lines[i].strip().split()[1])
  
    lines = lines[n_atoms-2:]
  
    for i in range(n_atoms-3):
      D[i] = float(lines[i].strip().split()[1])
  
    lines = lines[n_atoms-3:]
  modredundant = []
  i = 1
  line = lines[i].strip()
  while line[:4] != 'Grad':
    if len(line) > 0 and line.split()[0] == 'D':
      modredundant.append(line.split()[1:5])
    i = i+1
    line = lines[i].strip()

  modredundant = np.array(modredundant, dtype='int') - 1

  lines = lines[i:]
  lines = lines[9:]

  params_R = []
  params_A = []
  params_D = []
  params_value = []

  i = 0
  line = lines[i].strip()
  while line[0] == '!':
    line = line.strip('!').split()
    temp = line[1][1:].strip('()').split(',')
    if line[0][0] == 'R':
      params_R.append(np.array(temp, dtype='int'))
    elif line[0][0] == 'A':
      params_A.append(np.array(temp, dtype='int'))
    elif line[0][0] == 'D':
      params_D.append(np.array(temp, dtype='int'))
    params_value.append(float(line[2]))
    i = i + 1
    line = lines[i].strip()

  params_R = np.array(params_R) - 1
  params_A = np.array(params_A) - 1
  params_D = np.array(params_D) - 1
  params_value = np.array(params_value, dtype='float64')

  n_params = int(i)
  lines = lines[i+2:]

  n_opt = int(lines[0].strip().split()[-1])
  lines = lines[4:]

  for i in np.arange(n_opt)+1:
    lines = ParseOneOpt(i, lines)
