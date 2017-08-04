#!/usr/bin/env python
# coding: utf-8

from __future__ import division
import numpy as np
import scipy as sp
import sys,os

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
    for line in lines[1:1+noa]:
        e, x, y, z = line.split()
        elems.append(e)
        coords.append([x,y,z])
    return elems, np.array(coords, dtype=float)


#basis = load_basis('6-311g.basis')
basis = load_basis('3-21g.basis')
elems, coords = read_xyz('HeH.xyz')
#convert from A to Bohr
coords/=0.529177
cgf_list=[]
for i in range(len(elems)):
    e=elems[i]
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
        cgf_list.append((spdf,chemical_elems[elems[i]],new_coefs,coords[i])) 
        #need to end up summing over the primitives when there is a reference atom/radius r value
        # this should be done in a loop I think

                #Nuclear attraction matrix V
            #[[-6.01887445 -2.55025777 -0.35016867 -0.87003997]
            # [-2.55025777 -2.49407349 -0.87537411 -1.29768916]
            # [-0.35016867 -0.87537411 -2.77875362 -1.46151722]
            # [-0.87003997 -1.29768916 -1.46151722 -1.62944811]]

            #core with T+V as being evaluated
    
ref=0 # setting the atom name and coordinates of the respective atom - Rc
Zc=chemical_elems[elems[ref]]  
Rc=coords[ref]
# assuming all spdf = s
# finding number of basis functions
nbf = len(cgf_list)
# also want to consider the normalization factor
# either in put into integal or inclde them into the basis functions

## previous version for the overlap loops;
## wondering what is a more efficient way of looping over the basis wrt the future operations
overlap = np.zeros([nbf,nbf])
potential = np.zeros([nbf,nbf])
core = np.zeros([nbf,nbf])
kinetic = np.zeros([nbf,nbf])
for i in range(nbf):
    spdf1, Z1, coefs1, R1 = cgf_list[i]
    for j in range(i, nbf):
        spdf2, Z2, coefs2, R2 = cgf_list[j]
        S=0.0; V1=0.0; V2=0.0; T=0.0; H=0.0
        for alpha1,coef1 in coefs1:
            for alpha2,coef2 in coefs2:
                R=np.dot(R1-R2,R1-R2)
                val=(np.pi/(alpha1+alpha2))**(1.5)
                val*=np.exp((-alpha1*alpha2/(alpha1+alpha2))*R)
                val*=coef1*coef2  
                overlap[i,j]+=val
                ##now trying to calculate the kinetic energy
                ##Eq. A.11
                val=(alpha1*alpha2/(alpha1+alpha2))
                val*=(3.-2.*alpha1*alpha2/(alpha1+alpha2)*(R))
                val*=(np.pi/(alpha1+alpha2))**(1.5)
                val*=np.exp((-alpha1*alpha2/(alpha1+alpha2))*(R))
                val*=coef1*coef2  
                kinetic[i,j]+=val
                #**STRUGGLING HERE:
                #evaluating V uclear attraction integral A.33
                #need to have another loop to loop over all the Rc's and Zcs.
                Rp=(alpha1*R1+alpha2*R2)/(alpha1+alpha2)
                Rpc1=np.dot(Rp-R1,Rp-R1)
                if abs(Rpc1)<=1.0e-33:
                    value0=0
                else:
                    val=(-2*np.pi/(alpha1+alpha2))*Z1
                    val*=np.exp((-alpha1*alpha2/(alpha1+alpha2))*(R))
                    time=((alpha1+alpha2)*(Rpc1))
                    val*=(0.5*(np.pi/time)**(0.5))
                    val*=sp.special.erf(time**0.5) 
                    val*=coef1*coef2  
                    potential[i,j]+=val
                Rpc2=np.dot(Rp-R2,Rp-R2)
                if abs(Rpc2)<=1.0e-33:
                    value0=0
                else:
                    val=(-2*np.pi/(alpha1+alpha2))*Z2
                    val*=np.exp((-alpha1*alpha2/(alpha1+alpha2))*(R))
                    time=((alpha1+alpha2)*(Rpc2))
                    val*=(0.5*(np.pi/time)**(0.5))
                    val*=sp.special.erf(time**0.5) 
                    val*=coef1*coef2  
                    potential[j,i]+=val
    overlap[:,i]=overlap[i,:]
    kinetic[:,i]=kinetic[i,:]
    #potential[:,i]=potential[i,:]
    #core[:,i]=core[i,:]
print "Overlap matrix S"
print overlap
print "Kinetic energy matrix T"
print kinetic
print "Nuclear attraction matrix V"
print potential
print "Core Hamiltonian matrix H = T + V"
print core

#2 electron - 4 and 4 rather than 2 and 2 for loops


#
#Nuclear attraction matrix V
#[[-6.01887445 -2.55025777 -0.35016867 -0.87003997]
# [-2.55025777 -2.49407349 -0.87537411 -1.29768916]
# [-0.35016867 -0.87537411 -2.77875362 -1.46151722]
# [-0.87003997 -1.29768916 -1.46151722 -1.62944811]]
#
#Core Hamiltonian matrix H = T + V
#[[-2.10274492 -1.97130795 -0.42270654 -0.80566668]
# [-1.97130795 -1.91958399 -0.77842661 -1.14530694]
# [-0.42270654 -0.77842661 -1.22934868 -1.16836542]
# [-0.80566668 -1.14530694 -1.16836542 -1.35466011]]
#
#Two-electron Integrals G
#[[[[ 1.80325987  0.87650874  0.11829416  0.29968879]
#   [ 0.87650874  0.91700963  0.23781881  0.43889529]
#   [ 0.11829416  0.23781881  0.52774164  0.33692221]
#   [ 0.29968879  0.43889529  0.33692221  0.46613298]]
#
#  [[ 0.87650874  0.45207792  0.0628695   0.15616976]
#   [ 0.45207792  0.51467845  0.13770134  0.25089583]
#   [ 0.0628695   0.13770134  0.31290548  0.19886278]
#   [ 0.15616976  0.25089583  0.19886278  0.27375148]]
#
#Overlap matrix S
#[[ 1.00000146  0.59521639  0.09090182  0.21271836]
# [ 0.59521639  1.          0.33445798  0.58142644]
# [ 0.09090182  0.33445798  1.00000025  0.64589902]
# [ 0.21271836  0.58142644  0.64589902  1.        ]]
#
#Kinetic energy matrix T
#[[ 3.91612953  0.57894982 -0.07253787  0.06437328]
# [ 0.57894982  0.5744895   0.0969475   0.15238223]
# [-0.07253787  0.0969475   1.54940493  0.2931518 ]
# [ 0.06437328  0.15238223  0.2931518   0.274788  ]]

#  ...

#overlap = np.zeros([nbf,nbf])
#for i in range(nbf):
#    bas1=cgf_list[i][1]
#    for j in range(i, nbf):
#        bas2=cgf_list[j][1]
#        R=np.sum(np.square(cgf_list[i][2]-cgf_list[j][2]))
#        # going to compute 0,0, then 0,1, then 0,2...
#        # going to compute the upper triangle of the matrix
##        overlap[i,j] = Overlap_CGF(cgf_list[i],cgf_list[j])
#        for primrow in range(len(bas1[1])):
#            for primcol in range(len(bas2[1])):
#                overlap[i][j]+=bas1[primrow][1]*np.exp(-bas1[primcol][0]*R)
#                print overlap
#
#        # looking at equation A9 in the Szabo Ostlund book
#        # in this derivation, there is only alpha, but in final result, need to multiply by c's
#        # need to loop over primitive gaussians
#print overlap


##LOOK HERE FOR MORE INFO ON FUNCTIONS***


#Smatrix=np.zeros((nbf,nbf))
#Kmatrix=np.zeros((nbf,nbf))
##double loops
#for molrow in range(len(molecule)):
#    #columns
#    for molcol in range(len(molecule)):
#        #loop over the number of different primitives in the basis
#        #curious of a case of varying number of primitives used per matrix element?...I don't think so atm --changed my mind!/need to double loop over bases primitives and option of different number of primitives..
#        for primrow in range(len(STO3G)):
#            # need to have an empty matrix to start filling, that'd be the nxn
#            # now to start calling out the int function and then
#            # sum up the values in the matrix elements (I think?...)
#            for primcol in range(len(STO3G)):
#                Smatrix[molrow][molcol]+=intSAB(STO3G[primrow],STO3G[primcol],molecule[molrow][1],molecule[molcol][1])
#                Kmatrix[molrow][molcol]+=intKAB(STO3G[primrow],STO3G[primcol],molecule[molrow][1],molecule[molcol][1])
#print "Overlap Matrix Values:\n", Smatrix
#print "Kinetic Energy Matrix Values:\n", Kmatrix
#
