#!/usr/bin/env python
# coding: utf-8

from __future__ import division
import numpy as np
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
        cgf_list.append((spdf,new_coefs,coords[i])) 
        #need to end up summing over the primitives when there is a reference atom/radius r value
        # this should be done in a loop I think



def Overlap_CGF(bas1,bas2):
    result=0.0
    KEs=0.0
    spdf1, coefs1, R1 = bas1
    spdf2, coefs2, R2 = bas2
    for alpha1,coef1 in coefs1:
        for alpha2,coef2 in coefs2:
            R=np.dot(R1-R2,R1-R2)
            val=(np.pi/(alpha1+alpha2))**(1.5)
            val*=np.exp((-alpha1*alpha2/(alpha1+alpha2))*R)
            val*=coef1*coef2  
            result+=val
            ##now trying to calculate the kinetic energy
            #now evaluate the kinetic energy integral
            #da_a -- d and alpha values for atom A, positions 0 and 1
            #evaluating the overlap matrix from S&A:
            ##Eq. A.11
            #need normalization? this from last time? the following line is pretty sure as normalization
            val=(4.*alpha1*alpha2/(np.pi**2))**(0.75)
            val*=(np.pi/(alpha1+alpha2))**(1.5)
            val*=(alpha1*alpha2/(alpha1+alpha2))*(3.-2.*alpha1*alpha2/(alpha1+alpha2)*(R))
            val*=np.exp((-alpha1*alpha2/(alpha1+alpha2))*(R))
            KEs+=val
            #not need this? - from overlap matrix calculation
            #val*=da_a[0]*da_b[0]            
            # overlapof 2 
            # can even calculate overlap and kinetic energy
#        nbf=len(bastot)
#        for prim in range(nbf):
#            overlap+=bastot[prim][1]*np.exp(-bastot[prim][0]*R) #not quite accurate, need R variable based on location
        # this should actually be going over the overlap matrix such as what was in the previous code, where there are the factors in the front....        
    return KEs,result
    
    
# assuming all spdf = s
# finding number of basis functions
nbf = len(cgf_list)
# also want to consider the normalization factor
# either in put into integal or inclde them into the basis functions

## previous version for the overlap loops;
## wondering what is a more efficient way of looping over the basis wrt the future operations
overlap = np.zeros([nbf,nbf])
kinetic = np.zeros([nbf,nbf])
for i in range(nbf):
    for j in range(i, nbf):
        # going to compute 0,0, then 0,1, then 0,2...
        # going to compute the upper triangle of the matrix
        #not sure if want the two lists, or submit a total list of all the options; consider this later if it comes up
        # nice so far with only 2 atoms of interest to have the two entries for the two sets of basis functions; 
        # --would be nicer to be more general, could make this generalization in the CGF overlap function
        kinetic[i,j],overlap[i,j] = Overlap_CGF(cgf_list[i],cgf_list[j])
#        overlap[i,j] = Overlap_CGF(cgf_list[i],cgf_list[j])

        # looking at equation A9 in the Szabo Ostlund book
        # in this derivation, there is only alpha, but in final result, need to multiply by c's
        # need to loop over primitive gaussians
#    kinetic[:,i]=kinetic[i,:]
    overlap[:,i]=overlap[i,:]
print kinetic
print overlap

#2 electron - 4 and 4 rather than 2 and 2 for loops



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