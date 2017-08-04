# -*- coding: utf-8 -*-
"""
Created on Wed May 24 17:19:51 2017

@author: Lisa
"""
import numpy as np

# inputs
# 1) list of atoms and their positions (H, He)
# 2) An atomic basis set

# outputs
# 1) Build a hamiltonian matrix
# 2) diagonalize hamiltonian to obtain a set of 1-electron wavefunctions (orbitals) and their energies

# 3 elements
# each elemetn a list of 2 eleemtns
# one element exponent, the other a coefficient
#coeff, exp

#list of basis functions
#--centers
#--primitive gaussians
#--write loop over matrix
#-----loop over primitvie gaussians, get exponenets

#made a tuple, first with position, second of basis, within basis has the sets of 
# not really too sure how to make the library to call out the type of basis functions...

  
STO3G = [[0.444635,0.168856],[0.535328,0.623913],[0.154329,3.42525]]
molecule = [["H",[0,0,0]],["H",[1.4,0,0]]]
n=len(molecule)
Smatrix=np.zeros((n,n))
Kmatrix=np.zeros((n,n))

#da_a -- d and alpha values for atom A, positions 0 and 1
#evaluating the overlap matrix from S&A:
##Eq. A.9
def intSAB(da_a,da_b,Ra,Rb):
    Ra=np.array(Ra)
    Rb=np.array(Rb)
    # I don't remember what this first value means; maybe normalization?
    val=(4.*da_a[1]*da_b[1]/(np.pi**2))**(0.75)
    #the below are from S&A 
    val*=(np.pi/(da_a[1]+da_b[1]))**(1.5)
    #is this from A.7?
    val*=np.exp((-da_a[1]*da_b[1]/(da_a[1]+da_b[1]))*(np.sum(np.square(Ra-Rb))))
    val*=da_a[0]*da_b[0]
    return val

#now evaluate the kinetic energy integral
##Eq. A.11
def intKAB(da_a,da_b,Ra,Rb):
    Ra=np.array(Ra)
    Rb=np.array(Rb)
    #need normalization? this from last time?
    #val=(4.*da_a[1]*da_b[1]/(np.pi**2))**(0.75)
    val=(da_a[1]*da_b[1]/(da_a[1]+da_b[1]))*(3.-2.*da_a[1]*da_b[1]/(da_a[1]+da_b[1])*(np.sum(np.square(Ra-Rb))))
    val*=(np.pi/(da_a[1]+da_b[1]))**(1.5)
    val*=np.exp((-da_a[1]*da_b[1]/(da_a[1]+da_b[1]))*(np.sum(np.square(Ra-Rb))))
    #not need this? - from overlap matrix calculation
    #val*=da_a[0]*da_b[0]
    return val

    
#double loops
#for loops basis row, basis column:
#    for loops loop over 3 primitives of sto3g for each of the two atoms of interest, or however many atoms of interest
#rows 
for molrow in range(len(molecule)):
    #columns
    for molcol in range(len(molecule)):
        #loop over the number of different primitives in the basis
        #curious of a case of varying number of primitives used per matrix element?...I don't think so atm --changed my mind!/need to double loop over bases primitives and option of different number of primitives..
        for primrow in range(len(STO3G)):
            # need to have an empty matrix to start filling, that'd be the nxn
            # now to start calling out the int function and then
            # sum up the values in the matrix elements (I think?...)
            for primcol in range(len(STO3G)):
                Smatrix[molrow][molcol]+=intSAB(STO3G[primrow],STO3G[primcol],molecule[molrow][1],molecule[molcol][1])
                Kmatrix[molrow][molcol]+=intKAB(STO3G[primrow],STO3G[primcol],molecule[molrow][1],molecule[molcol][1])
print "Overlap Matrix Values:\n", Smatrix
print "Kinetic Energy Matrix Values:\n", Kmatrix





