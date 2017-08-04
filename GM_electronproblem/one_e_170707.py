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

def Overlap_CGF(alpha1,alpha2,coefs,R):            
    val=(np.pi/(alpha1+alpha2))**(1.5)
    val*=np.exp((-alpha1*alpha2/(alpha1+alpha2))*R)
    val*=coefs  
    return val
            
def Kinetic_CGF(alpha1,alpha2,coefs,R):
    ##Eq. A.11
    val=(alpha1*alpha2/(alpha1+alpha2))
    val*=(3.-2.*alpha1*alpha2/(alpha1+alpha2)*(R))
    val*=(np.pi/(alpha1+alpha2))**(1.5)
    val*=np.exp((-alpha1*alpha2/(alpha1+alpha2))*(R))
    val*=coefs  
    return val

#def Potential_CGF(alpha1,alpha2,coef1,coef2,R,Rc,Zc):
#    #evaluating V uclear attraction integral A.33  
#    val=(-2*np.pi/(alpha1+alpha2))*Zc
#    val*=np.exp((-alpha1*alpha2/(alpha1+alpha2))*(R))
#    val*=(0.5*(np.pi/time)**(0.5))
#    val*=sp.special.erf(time**0.5) 
#    val*=coefs  
#    return val
    
# assuming all spdf = s
nbf = len(cgf_list)
overlap = np.zeros([nbf,nbf])
kinetic = np.zeros([nbf,nbf])
potential = np.zeros([nbf,nbf])
for i in range(nbf):
    spdf1, Z1, coefs1, R1 = cgf_list[i]
    for j in range(i,nbf):
        spdf2, Z2, coefs2, R2 = cgf_list[j]
        R=np.dot(R1-R2,R1-R2)
        for alpha1,coef1 in coefs1:
            for alpha2,coef2 in coefs2:
                coefs=coef1*coef2
                S=Overlap_CGF(alpha1,alpha2,coefs,R)
                overlap[i,j]+=S
                T=Kinetic_CGF(alpha1,alpha2,coefs,R)
                kinetic[i,j]+=T
                Rp=(alpha1*R1+alpha2*R2)/(alpha1+alpha2)
                #need to have another loop to loop over all the Rc's and Zcs.
                for k in range(len(coords)):
                    Zc=chemical_elems[elems[k]]
                    Rpc=np.dot(Rp-coords[k],Rp-coords[k])
                    if abs(Rpc)==0.0:
                        Fo=1
                    else:
                        time=((alpha1+alpha2)*(Rpc))
                        Fo=(0.5*(np.pi/time)**(0.5))*sp.special.erf(time**0.5)                        
                    V=(-2*np.pi/(alpha1+alpha2))*Zc
                    V*=Fo
                    V*=np.exp((-alpha1*alpha2/(alpha1+alpha2))*(R))
                    V*=coefs
                    potential[i,j]+=V
    overlap[:,i]=overlap[i,:]
    kinetic[:,i]=kinetic[i,:]
    potential[:,i]=potential[i,:]
    core=kinetic+potential
print "Overlap matrix S"
print overlap
print "Kinetic energy matrix T"
print kinetic
print "Nuclear attraction matrix V"
print potential
print "Core Hamiltonian matrix H = T + V:"
print core

S=overlap
H=core
T=kinetic
V=potential

#eigh where h is for hermitian
Eval,Evec=np.linalg.eigh(S)
s=np.diag(Eval**-0.5)
U=np.zeros(np.shape(Evec))
for i in range(len(Evec)):
    U[:,i]=Evec[i]
#[np.argsort(Eval)]
val=np.dot(s,U.T)
S_12=np.dot(U,val)
val=np.dot(H,S_12)
Hp=np.dot(S_12,val)
e,Cp=np.linalg.eig(Hp)
C=np.dot(S_12,Cp)
print "Coefficients:"
print C

#2 electron - 4 and 4 rather than 2 and 2 for loops

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
