# rigorous_conley_index.py

import traceback
import ctypes
import numpy as np
from collections import OrderedDict
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os.path

import matplotlib
#matplotlib.use('Agg') 
import matplotlib.pyplot as plt

import cellComplex_python as cellComplex


class ConleyExample(object):
    """ A ConleyExample consists of a VectorField and a
        CandidateIsolatingBlock.  It has the capability of validating
        candidate isolating blocks and computing the homological
        conley index of a validated isolating block. """
    def __init__(self,example_name):
        self.example_name = example_name
        self.vector_field = VectorField(self.example_name)
        self.candidate_block = None

    def set_candidate_block(self,candidate_block_name):
        pass

    def validate_candidate_block(self):
        exit_flag = 
        cell_complex = 
        return exit_flag,candidate_block
    
class VectorField(object):
    """ A vector field has a callable expression, equilibria,
        jacobians, and various default settings for plotting and
        integrating. """
    def __init__(self,example_name):
        self.example_name = example_name
    
class CandidateIsolatingBlock(object):
    """ A candidate isolating block has a parameterization and can be
        rigorously validated and if valid, one can compute the
        homological conley index of the block. """
    def __init__(self,block_name,parameterization):
        self.block_name = block_name
        self.parameterization = parameterization
        self.is_validated = False
       
    def set_parameterization(self,block_name):
        parameterization = lib.get_candidate_block_parameterization(self.example_name, candidate_block_name)
   
    def rigorously_validate_block(self):
       pass
    
    def compute_homological_conley_index(self,candidate_block):
        pass












def getUserDefinedFunction():
    #{{{
    return lib.getUserDefinedFunction()
    #}}}

def get_exit_set_function(example,block):
    #{{{
    return lib.get_exit_set_function()
    #}}}

def get_tangency_test_function():
    #{{{
    return lib.get_tangency_test_function()
    #}}}

def getPlottableData2d(p,Nx=200,Ny=200):
    #{{{
    a,b = p[0]
    c,d = p[1]
    hx = float(b-a)/(Nx-1)
    hy = float(d-c)/(Ny-1)
    z = np.zeros((Nx,Ny))

    for i in range(Nx):
        x = a+i*hx
        for j in range(Ny):
            y = c+j*hy
            fVal = lib.pointEvalUserDefinedFunction_2d(x,y)
            #print fVal
            z[i,j] = fVal
    xx,yy = np.meshgrid(np.linspace(a,b,Nx),np.linspace(c,d,Ny))
    return (xx,yy,z)
    #}}}

def getPlottableData3d(x,y,z):
    #{{{
    fVal = lib.pointEvalUserDefinedFunction_3d(x,y,z)
    return fVal
    #}}}

