#!/usr/bin/env python
# coding: utf-8

from __future__ import division
import numpy as np

chemical_elems = {'H':1, 'He':2}
chemical_elems['Li'] = 3

def load_basis(filename):
    basis = dict()
    reading = False
    for line in open(filename):
        if line.startswith('#'):
            continue
        elif line.startswith('BASIS'):
            reading = True
        elif line.startswith('END'):
            break
        elif reading is True:
            ls = line.split()
            if len(ls) == 2 and ls[0].isalpha():
                elem, spdf = ls
                if elem not in basis:
                    basis[elem] = []
                coefs = list()
                basis[elem].append((spdf, coefs))
            else:
                v = [float(i) for i in ls]
                coefs.append(v)
    return basis
        
def read_xyz(filename):
    lines = open(filename).readlines()
    noa = int(lines[0])
    elems = []
    coords = []
    for line in lines[2:2+noa]:
        e, x, y, z = line.split()
        elems.append(e)
        coords.append([x,y,z])
    return elems, np.array(coords, dtype=float)


basis = load_basis('6-311g.basis')
elems, coords = read_xyz('HeH.xyz')
Ra = coords[0]
# this part of the code didnt loop over the elements, but only took the first element
cgf_list = []
##NEED TO CONSIDER WHAT TO INCLUDE IN THIS LIST****
e = elems[0]
#H is an example return value
basis_info = basis[e]
#return value: orbital, [(alpha and c values for the cartesian gaussian basis functions)]
#[('S', [[33.865, 0.0254938], [5.09479, 0.190373], [1.15879, 0.852161]]), ('S', [[0.32584, 1.0]]), ('S', [[0.102741, 1.0]])]
for spdf, coefs in basis_info:
    new_coefs = []
    for alpha, c in coefs:
        norm = (2.*alpha/np.pi)**0.75
        c*=norm
        new_coefs.append([alpha,c])
    # include the c*norm as well as the PGF values (not including the norm'd, or var's needed)
    #orb, alpha, norm'd c
    cgf_list.append((spdf,new_coefs,Ra))
    
#first one spdf, next one is coefficients
# spherical only (5?) lz only 5 differnt compoents
# in carteisan fraome 10 different f
# when listing all the things that can add up to 3
# only 6 f functions that can add to 3
# since only working wisth s fucnations, not worry about expanding the same coeffs into px py pz
# for higher angular momentums, need to expand into higher step

# going to assume all spdf = s
nbf = len(cgf_list)
# also want to consider the normalization factor
# either in put into integal or inclde them into the basis functions
print cgf_list

overlap = np.zeros([nbf,nbf])
for i in range(nbf):
    for j in range(i, nbf):
        # going to compute 0,0, then 0,1, then 0,2...
        # going to compute the upper triangle of the matrix
        ovelap[i,j] = Overlap_CGF(cgf_list[i],cgf_list[j])
        # looking at equation A9 in the Szabo Ostlund book
        # in this derivation, there is only alpha, but in final result, need to multiply by c's
        # need to loop over primitive gaussians

#def Overlap_CGF(cgf_list[i],cgf_list[j]):
#    result = 0
#    spdf1, coefs1, Ra1 = cgf_list[i]
#    # want to make sure that the the constant values are not labeled the same such that overwriting
#    spdf2, coefs2, Ra2 = cgf_list[j]
#    for pgf1:
#        coefs1
#        for pgf2:
#            coefs2
#            # now can easily insert the formulas
#            # also need a distance between the two; there's a hint that you can directly do Ra - Rb
#            np.dot(Ra-Rb,Ra-Rb)
#            # overlapof 2 
## can even calculate overlap and kinetic energy
