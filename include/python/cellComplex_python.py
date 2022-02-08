# cellComplex_python.py

import traceback
import ctypes
import numpy as np
from collections import OrderedDict
from mpl_toolkits.mplot3d.art3d import PolyCollection, Poly3DCollection
import os.path

import matplotlib
#matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# what a cludge...
import inspect
class DummyClass: pass
dirname = os.path.dirname(os.path.abspath(inspect.getsourcefile(DummyClass)))
print('looking in: ', dirname)
lib = ctypes.cdll.LoadLibrary('%s/cellComplex_lib.so' %(dirname))
# end cludge

__all__=['CellComplex']

POINT = [[-1,1],[-1,1],[-1,1]] # for plotting later...
del_color = 'cyan'    ## negative
undel_color = 'blue'  ## positive
mixed_boundary_color = 'red'  ## 
#mixed_boundary_color = del_color  ## 

def get_userDefinedFunction():
    #{{{
    return lib.get_userDefinedFunction()
    #}}}

def get_exit_set_function():
    #{{{
    return lib.get_exit_set_function()
    #}}}

def get_tangency_test_function():
    #{{{
    return lib.get_tangency_test_function()
    #}}}

# def get_energy_density_function():
#     #{{{
#     return lib.get_energy_density_function()
#     #}}}


def get_block_parameterization(parameterization='A'):
    #{{{
    return lib.get_block_parameterization()
    #}}}

def get_vector_field():
    #{{{
    return lib.get_vector_field()
    #}}}

def eval_autograd(p=None,func=None):

    autograd = -9 # sentinal value, autograd should be a vector of intervals
    pair = ctypes.c_double * 2
    point = (pair * len(p))(*(pair(*q[:2]) for q in p))

    if func == 'exit_set':
        f = get_exit_set_function()
    if func == 'tangency_test':
        f = get_tangency_test_function()
    if func == 'energy_density':
        f = get_energy_density_function()
    
    lib.eval_autograd.argtypes = \
        [(ctypes.c_double * 2)*len(p),ctypes.c_int,ctypes.c_void_p]
    
    autograd = lib.eval_autograd(point,len(p),f)
    return autograd.contents

def eval_vector_field(x,y,z):
    #{{{
    fVal = lib.eval_vector_field(x,y,z)
    return fVal
    #}}}

def eval_exit_set(p,Nx=200,Ny=200):
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
            fVal = lib.eval_exit_set(x,y)
            #print fVal
            z[i,j] = fVal
    xx,yy = np.meshgrid(np.linspace(a,b,Nx),np.linspace(c,d,Ny))
    return (xx,yy,z)
    #}}}

def eval_tangency_test_w(p,Nx=200,Ny=200):
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
            fVal = lib.eval_tangency_test_w(x,y)
            #print fVal
            z[i,j] = fVal
    xx,yy = np.meshgrid(np.linspace(a,b,Nx),np.linspace(c,d,Ny))
    return (xx,yy,z)
    #}}}

def eval_tangency_test_v(p,Nx=200,Ny=200):
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
            fVal = lib.eval_tangency_test_v(x,y)
            #print fVal
            z[i,j] = fVal
    xx,yy = np.meshgrid(np.linspace(a,b,Nx),np.linspace(c,d,Ny))
    return (xx,yy,z)
    #}}}

def eval_func_on_surface(p=None,parameterization='A',func=None,Nx=200,Ny=200):
    #{{{
    a,b = p[0]
    c,d = p[1]

    xx  = np.linspace(a,b,Nx)
    yy  = np.linspace(c,d,Ny)

    txx = np.zeros((Nx,Ny))   
    tyy = np.zeros((Nx,Ny))   
    tzz = np.zeros((Nx,Ny))   

    fVal   = np.zeros((Nx,Ny))
    
    if func == 'exit_set':
        f_handle = lib.eval_exit_set
    if func == 'tangency_test':
        f_handle = lib.eval_tangency_test
    #if func == 'energy_density':
    #    f_handle = lib.eval_energy_density
    
    for i,x in enumerate(xx):
        for j,y in enumerate(yy):
            f = f_handle(x,y)  
            #f = lib.eval_exit_set(x,y)
            fVal[i,j] = f
            TX = lib.eval_block_parameterization(x,y)
            tx,ty,tz  = TX.contents

            txx[i,j],tyy[i,j],tzz[i,j] = tx,ty,tz  


    return (txx,tyy,tzz,fVal)
    #}}}

#def eval_block_parameterization(p,Nx=200,Ny=200):
#    #{{{
#    a,b = p[0]
#    c,d = p[1]
#    hx = float(b-a)/(Nx-1)
#    hy = float(d-c)/(Ny-1)
#    z = np.zeros((Nx,Ny))
#
#    for i in range(Nx):
#        x = a+i*hx
#        for j in range(Ny):
#            y = c+j*hy
#            fVal = lib.eval_block_parameterization(x,y)
#            #print fVal
#            z[i,j] = fVal
#    xx,yy = np.meshgrid(np.linspace(a,b,Nx),np.linspace(c,d,Ny))
#    return (xx,yy,z)
#    #}}}

#def getScalarData(p,func=None,Nx=200,Ny=200):
#    #{{{
#    a,b = p[0]
#    c,d = p[1]
#    hx = float(b-a)/(Nx-1)
#    hy = float(d-c)/(Ny-1)
#    z = np.zeros((Nx,Ny))
#
#    if func == 'userDefinedFunction':
#        f_handle = lib.pointEvalUserDefinedFunction
#    #if func == 'exit_set_u'
#    #    f_handle = lib.pointEval
#
#    for i in range(Nx):
#        x = a+i*hx
#        for j in range(Ny):
#            y = c+j*hy
#            fVal = f_handle(x,y)
#            #fVal = lib.pointEvalUserDefinedFunction_2d(x,y)
#            #print fVal
#            z[i,j] = fVal
#    xx,yy = np.meshgrid(np.linspace(a,b,Nx),np.linspace(c,d,Ny))
#    return (xx,yy,z)
#    #}}}

def getPlottableData2d(p,func=None,Nx=200,Ny=200):
    #{{{
    a,b = p[0]
    c,d = p[1]
    hx = float(b-a)/(Nx-1)
    hy = float(d-c)/(Ny-1)
    z = np.zeros((Nx,Ny))

    if func == 'exit_set':
        f_handle = lib.eval_exit_set
    elif func == 'tangency_test':
        f_handle = lib.eval_tangency_test
    #elif func == 'energy_density':
    #    f_handle = lib.eval_energy_density
    elif func == 'userDefinedFunction2d':
        f_handle = lib.pointEvalUserDefinedFunction_2d

    for i in range(Nx):
        x = a+i*hx
        for j in range(Ny):
            y = c+j*hy
            fVal = f_handle(x,y)
            #fVal = lib.pointEvalUserDefinedFunction_2d(x,y)
            #print fVal
            z[i,j] = fVal
    xx,yy = np.meshgrid(np.linspace(a,b,Nx),np.linspace(c,d,Ny))
    return (xx,yy,np.transpose(z))
    #}}}

def getPlottableData3d(x,y,z):
    #{{{
    #a,b = p[0]
    #c,d = p[1]
    #e,f = p[2]
    #hx = float(b-a)/(Nx-1)
    #hy = float(d-c)/(Ny-1)
    #hz = float(e-f)/(Nz-1)
    #xx = np.zeros((Nx,Ny,Nz))
    #yy = np.zeros((Nx,Ny,Nz))
    #zz = np.zeros((Nx,Ny,Nz))
    #w  = np.zeros((Nx,Ny,Nz))

    #for i in range(Nx):
    #    x = a+i*hx
    #    for j in range(Ny):
    #        y = c+j*hy
    #        for k in range(Nz):
    #            z = e+k*hz
    #            xx[i,j,k] = x
    #            yy[i,j,k] = y
    #            zz[i,j,k] = z
    #            fVal = lib.pointEvalUserDefinedFunction_3d(x,y,z)
    #            #print fVal
    #            w[i,j,k] = fVal
    ##xx,yy,zz = np.meshgrid(np.linspace(a,b,Nx),np.linspace(c,d,Ny),np.linspace(e,f,Nz))
    ##xx, yy, zz = np.mgrid[a:b:Nx, c:d:Ny, e:f:Nz]
    
    fVal = lib.pointEvalUserDefinedFunction_3d(x,y,z)
    return fVal

    #}}}

def validateDomainByValidatingTopDimensionalCells(p=None,f=None):
    #{{{
    #global POINT 
    POINT = p # for plotting later...
    pair = ctypes.c_double * 2
    point = (pair * len(p))(*(pair(*q[:2]) for q in p))
    lib.validateDomainByValidatingTopDimensionalCells.argtypes = \
        [(ctypes.c_double * 2)*len(p),ctypes.c_int,ctypes.c_void_p]
    validatedCmplx = lib.validateDomainByValidatingTopDimensionalCells(point,len(p),f)

    return validatedCmplx
    #}}}

def is_positive_on_entire_cell(p,f):
    #{{{
    pair = ctypes.c_double * 2
    point = (pair * len(p))(*(pair(*q[:2]) for q in p))
    lib.is_positive_on_entire_cell.argtypes = \
        [(ctypes.c_double * 2)*len(p),ctypes.c_int,ctypes.c_void_p]
    is_positive = lib.is_positive_on_entire_cell(point,len(p),f)

    return is_positive
    #}}}

def computeNodalDomainOfAFunction(p,f):
    #{{{
    global POINT 
    POINT = p # for plotting later...
    pair = ctypes.c_double * 2
    point = (pair * len(p))(*(pair(*q[:2]) for q in p))
    lib.computeNodalDomain.argtypes = \
        [(ctypes.c_double * 2)*len(p),ctypes.c_int,ctypes.c_void_p]
    validatedCmplx = lib.computeNodalDomain(point,len(p),f)
    return validatedCmplx
    #}}}

def computePersistenceOfFunction(p):
    #{{{
    global POINT 
    POINT= p # for plotting later...
    pair = ctypes.c_double * 2
    point = (pair * len(p))(*(pair(*q[:2]) for q in p))
    lib.computePersistence.argtypes = [(ctypes.c_double * 2)*len(p),ctypes.c_int]
    validatedCmplx_withPersist = lib.computePersistence(point,len(p))
    return validatedCmplx_withPersist
    #}}}

def isoblockval(p=None,exit_set=None,tangency_test=None):
    #{{{
    """ deleted cells are negative, undeleted cells are positive """ 
    
    is_valid_block = False # output of this boolean valued function
    max_depth   = 10    # how many times to subdivide individual cell

    # get CW complex approximation to nodal domains of the exit set
    Q_orig = validateDomainByValidatingTopDimensionalCells(p=p,f=exit_set)
    print('exit_set cw_complex validated, now going for isolating block validation')
    Q = Q_orig.top_dimensional_elemen
    T,V = [],[]
    
    while Q:
    
        for q in Q:
            if q.has_mixed_boundary:
                T.append(q)
            else:
                V.append(q)
        Q = []

        while T:
            t = T.pop() ## t is a top dimensional cell having mixed boundary from Q: nodal domains of exit_set_u 
            
            if True :# is_positive_on_entire_cell(t.coords,tangency_test):
                V.append(t)
            else:
                cmplx_t = t.make_cellComplex()
                cell_at_top = cmplx_t.top_dimensional_elemen[0]
                
                dim_of_cell_at_top = cell_at_top.dim
                cell_struct = cmplx_t.divideCubeAlongLongestAxis(cmplx_t.top_dimensional_elemen[0])
                t1 = cell_struct.cell_one
                t2 = cell_struct.cell_two
                A = validateDomainByValidatingTopDimensionalCells(p=t1.coords,f=exit_set)
                B = validateDomainByValidatingTopDimensionalCells(p=t2.coords,f=exit_set)
                A = A.top_dimensional_elemen
                B = B.top_dimensional_elemen
                Q.extend(A)   # list.extend() appends multiple items to list at once
                Q.extend(B)

    if not T: ## w is positive on every cell having mixed boundary in nodal domains of exit_set_u
        is_valid_block = True

    print('block is a validated isolating block, handing back to caller')
    return Q_orig,V,is_valid_block

    """ NOTES: 
                This is essentially the algorithm from the Wanner
                Stephens SIADS 2014 paper.  Inside the "while T:" loop, we want to
                divide a Cell t \in T if tangency_test is not positive on t, however
                the original divideCube_...() functions in the C++ nodal domain
                validation data structures are defined for CellComplex data types only, not
                for Cell data types. Hence, we create the seemingly extraneous  
                CellComplex cmplx_t from the Cell t we are interested in, then
                divide the complex, and then extract the Cells t1 and t2.
        TODO:   
                Aug 12, 2015: Implement a maximum subdivision depth for Cells which
                had has_mixed_boundary==True and for which
                tangency_test is not positive_on_entire_cell.  
    """
    #}}}


class Cell(ctypes.c_void_p):
   
    def __init__(self,p):
        #{{{
        self.p = p
        pair = ctypes.c_double * 2
        point = (pair * len(self.p))(*(pair(*q[:2]) for q in self.p))
        lib.new_cell.argtypes = [(ctypes.c_double * 2)*len(self.p), ctypes.c_int]
        self.value = lib.new_cell(point, len(self.p)).value
        #}}}
    
    @property
    def dim(self):
        #{{{
        return lib.dim_cell(self)
        #}}}
    
    @property
    def coords(self):
        #{{{ 
        dimOfContainingComplex = lib.get_dimOfContainingComplex(self) # ACK!!
        coord = []
        for i in range(dimOfContainingComplex):
            c = lib.get_coordsOfCellAtLoc(self,i)
            coord.append([c.contents[0],c.contents[1]])
        return coord
        #}}}
    
    def plot(self,ax,showDelCells=[0,1,2,3],showUndelCells=[0,1,2,3],
             showMixedCells=[0,1,2,3],
             total_dimension=2, mixed_bd_color=False, 
             undetermined_evaluation=False, c='k',ms=3,lw=5,
             pos_face_alpha=0.15,neg_face_alpha=0.15,
             mixed_bd_face_alpha=0.15,
             **kwargs):                  
        #{{{  
        isDeleted = lib.isDeleted(self)
        has_mixed_boundary = False
        coords = self.coords
        verts = sorted(get_vertsOfCell(self))
        verts = [v[0:total_dimension] for v in verts] 

        
        dim = self.dim
        if isDeleted:
            #print "is deleted"
            kwargs_vertex = {'color':del_color,'marker':'o','mec':del_color,
                             'markersize':ms}
            kwargs_edge   = {'color':del_color,'lw':lw}
            #kwargs_face   = {'facecolor':'cyan','alpha':0.15} # , 'edgecolor':'cyan'}
            kwargs_face   = {'facecolor':del_color,'alpha':neg_face_alpha ,'edgecolor':del_color}
            kwargs_cube   = {'facecolor':del_color,'alpha':1.0} # , 'edgecolor':'None'}
        else:
            #print "is not deleted"
            kwargs_vertex = {'color':undel_color,'marker':'o','mec':undel_color,
                             'markersize':ms}
            kwargs_edge   = {'color':undel_color,'lw':lw}
            #kwargs_face   = {'facecolor':'blue','alpha':0.5} #,'edgecolor':'blue'}
            kwargs_face   = {'facecolor':undel_color,'alpha':pos_face_alpha,'edgecolor':undel_color}
            kwargs_cube   = {'facecolor':undel_color,'alpha':1.0} #,'edgecolor':'None'}
        
        if self.has_mixed_boundary:
            has_mixed_boundary=True
            if mixed_bd_color:
                mixed_boundary_color = mixed_bd_color
            #print "has mixed boundary"
            kwargs_vertex = {'color':mixed_boundary_color,'marker':'o','mec':'k',
                             'markersize':ms}
            kwargs_edge   = {'alpha':1.0,'color':mixed_boundary_color,'lw':lw}
            #kwargs_face   = {'facecolor':'blue','alpha':0.5} #,'edgecolor':'blue'}
            kwargs_face   = {'facecolor':mixed_boundary_color,'alpha':mixed_bd_face_alpha,'edgecolor':undel_color}
            kwargs_cube   = {'facecolor':mixed_boundary_color,'alpha':1.0} #,'edgecolor':'None'}
        if undetermined_evaluation:
            #print ""
            kwargs_vertex = {'color':c,'alpha':0.2,'marker':'o','mec':'k',
                             'markersize':ms}
            kwargs_edge   = {'color':c,'alpha':0.2,'lw':lw}
            #kwargs_face   = {'facecolor':'blue','alpha':0.5} #,'edgecolor':'blue'}
            kwargs_face   = {'facecolor':c,'alpha':0.2,'edgecolor':'k'}
            kwargs_cube   = {'facecolor':c,'alpha':0.2} #,'edgecolor':'None'}

        
        if dim == 0:

            if isDeleted and dim in showDelCells:
                ax.plot(*[[v] for v in verts[0][0:total_dimension]],**kwargs_vertex)
            if not isDeleted and dim in showUndelCells:
                ax.plot(*[[v] for v in verts[0][0:total_dimension]],**kwargs_vertex)
            
            #except TypeError:
            #    ax.plot([verts[0][0]],[verts[0][1]],[verts[0][2]],**kwargs_vertex)

        elif dim == 3:
            f1 = Poly3DCollection(list(zip([verts[0]],[verts[2]],               
                                      [verts[3]],[verts[1]])),**kwargs_cube) 
            f2 = Poly3DCollection(list(zip([verts[2]],[verts[6]],               
                                      [verts[7]],[verts[3]])),**kwargs_cube) 
            f3 = Poly3DCollection(list(zip([verts[6]],[verts[4]],               
                                      [verts[5]],[verts[7]])),**kwargs_cube) 
            f4 = Poly3DCollection(list(zip([verts[4]],[verts[0]],               
                                      [verts[1]],[verts[5]])),**kwargs_cube) 
            f5 = Poly3DCollection(list(zip([verts[1]],[verts[3]],               
                                      [verts[7]],[verts[5]])),**kwargs_cube) 
            f6 = Poly3DCollection(list(zip([verts[0]],[verts[2]],               
                                      [verts[6]],[verts[4]])),**kwargs_cube) 
            
            if has_mixed_boundary and dim in showMixedCells:
                ax.add_collection3d(f1)
                ax.add_collection3d(f2)
                ax.add_collection3d(f3)
                ax.add_collection3d(f4)
                ax.add_collection3d(f5)
                ax.add_collection3d(f6)
            if isDeleted and dim in showDelCells:
                ax.add_collection3d(f1)
                ax.add_collection3d(f2)
                ax.add_collection3d(f3)
                ax.add_collection3d(f4)
                ax.add_collection3d(f5)
                ax.add_collection3d(f6)
            if not isDeleted and dim in showUndelCells:
                ax.add_collection3d(f1)
                ax.add_collection3d(f2)
                ax.add_collection3d(f3)
                ax.add_collection3d(f4)
                ax.add_collection3d(f5)
                ax.add_collection3d(f6)
        
        elif dim == 2:
            if total_dimension == 3:
                face = Poly3DCollection(list(zip([verts[0]],[verts[2]],               
                                            [verts[3]],[verts[1]])),**kwargs_face) 
                if has_mixed_boundary and dim in showMixedCells:
                    coll0 = ax.add_collection3d(face)

                if isDeleted and dim in showDelCells:
                    coll1 = ax.add_collection3d(face)
                    #if kwargs_face['facecolor']=='orange':
                    #    coll1.set_alpha(0.65)
                if not isDeleted and dim in showUndelCells:
                    coll2 =ax.add_collection3d(face)
                    #if kwargs_face['facecolor']=='orange':
                    #    coll2.set_alpha(0.65)
                    
            else:
                face = PolyCollection(list(zip([verts[0]],[verts[2]],               
                                            [verts[3]],[verts[1]])),**kwargs_face) 
                
                if has_mixed_boundary and dim in showMixedCells:
                    coll0 = ax.add_collection(face)
                if isDeleted and dim in showDelCells:
                    coll1 = ax.add_collection(face)
                    #if kwargs_face['facecolor']=='orange':
                    #    coll1.set_alpha(0.65)

                if not isDeleted and dim in showUndelCells:
                    coll2 = ax.add_collection(face)
                    #if kwargs_face['facecolor']=='orange':
                    #    coll2.set_alpha(0.65)


        elif dim == 1:
            if total_dimension == 3:
                if has_mixed_boundary and dim in showMixedCells:
                    ax.add_collection3d(Poly3DCollection(
                                        list(zip([verts[0]],[verts[1]])),**kwargs_edge))
                
                if isDeleted and dim in showDelCells:
                    ax.add_collection3d(Poly3DCollection(
                                        list(zip([verts[0]],[verts[1]])),**kwargs_edge))
                if not isDeleted and dim in showUndelCells:
                    ax.add_collection3d(Poly3DCollection(
                                        list(zip([verts[0]],[verts[1]])),**kwargs_edge))

            else:
                if has_mixed_boundary and dim in showMixedCells:
                    coll = ax.add_collection(PolyCollection(
                                             list(zip([verts[0]],[verts[1]])),**kwargs_edge))
                if isDeleted and dim in showDelCells:
                    coll = ax.add_collection(PolyCollection(
                                             list(zip([verts[0]],[verts[1]])),**kwargs_edge))
                    #if kwargs_edge['color'] == 'orange':
                    #    coll.set_alpha(0.9)
                if not isDeleted and dim in showUndelCells:
                    ax.add_collection(PolyCollection(
                                      list(zip([verts[0]],[verts[1]])),**kwargs_edge))

        
        #}}}

    def delete(self):
        #{{{
        self = lib.delet(self)
        #}}}

    def undelete(self):
        #{{{
        self = lib.undelet(self)
        #}}}

    @property
    def isDeleted(self,c):
        #{{{
        return lib.isDeleted(c)
        #}}}

    @property
    def has_mixed_boundary(self):
    #{{{
        has_mx_bd            = False
        contains_deleted     = False
        contains_undeleted   = False
         
        whole_boundary_of_cube = self.computeWholeBoundaryOfCube()

        while whole_boundary_of_cube:
            boundary_cell = whole_boundary_of_cube.pop()
            
            if lib.isDeleted(boundary_cell):
                contains_deleted = True
            if not lib.isDeleted(boundary_cell):
                contains_undeleted = True
            
            if contains_deleted and contains_undeleted:
                has_mx_bd = True
                break
        return has_mx_bd
    #}}}

    def make_cellComplex(self):
    #{{{
        # reconstruct p, i.e. invert coords()
        p = self.coords
        return CellComplex(p=p)
    #}}}

    def computeWholeBoundaryOfCube(self):
    #{{{
        whole_bd = []
        size = lib.size_wholeBoundaryOfCube(self)
        for j in range(size):
            whole_bd.append(lib.get_boundaryElementAtPosition(self,j))
        return whole_bd
    #}}}

class CellComplex(ctypes.c_void_p):
    
    def __init__(self,p=[[0,0]]):
        #{{{
        #while len(p) != 3:  
        #    p.append([0,0])
        pair = ctypes.c_double * 2
        point = (pair * len(p))(*(pair(*q[:2]) for q in p))
        self.value = lib.new_cellComplex(point, len(p)).value
        #}}}
        
    @property
    def dim(self):
        #{{{
        return lib.dim_cellComplex(self)
        #}}}

    @property
    def top_dimensional_elemen(self):
    #{{{
        size  = lib.size_elementsAtSpecifiedDim(self,self.dim-1)
        cells = []
        for j in range(size):
            cells.append(lib.get_elementAtSpecifiedDimAndLoc(self,self.dim-1,j))
        return cells
    #}}}

    @property
    def elemen(self):
        #{{{
        el = []
        for i in range(self.dim):
            size  = lib.size_elementsAtSpecifiedDim(self,i)
            cells = []
            for j in range(size):
                cells.append(lib.get_elementAtSpecifiedDimAndLoc(self,i,j))
            el.append(cells)
        return el
        #}}}
    
    @property
    def is_entirely_positive(self):
    #{{{
        is_entirely_pos = True
        
        print('about to call computeWholeBoundary from within python')

        whole_boundary_of_cube = self.computeWholeBoundaryOfCube()

        while whole_boundary_of_cube:
            boundary_cell = whole_boundary_of_cube.pop()
            
            if lib.isDeleted(boundary_cell):
                is_entirely_pos = False
                break
        
        return is_entirely_pos
    #}}}

    def divideTheCubeByGivenFactor(self,cube,n):
        #{{{ 
        self = lib.divideTheCubeByGivenFactor(self,cube,n)
        #}}}
    
    def divideCubeAlongLongestAxis(self,cell):
        #{{{ 
        two_cell_struct = lib.divideCubeAlongLongestAxis(self,cell)
        return two_cell_struct
        #}}}

    def divideCube(self,cell,n):
        #{{{ 
        two_cell_struct = lib.divideCube(self,cell,n)
        return two_cell_struct
        #}}}

    def plot(self,ax,legend=True,**kwargs):
        #{{{ 
        
        for i,cells in enumerate(self.elemen):
            print("number of cells of dim %d: %d" %(i,len(cells)))
            for c in cells:
                c.plot(ax,**kwargs)
        
        if legend:
            plot_legend(ax=ax,legend_u=True)
        
        #}}}
   
    def writeToFile(self,fname='../../data/testing.txt'):
        #{{{ 

        f = open(fname,'w')

        for i,cells in enumerate(self.elemen):
            print("number of cells of dim %d: %d" %(i,len(cells)))
            f.write('dim %d\n' %(i))
            for k,c in enumerate(cells):
                flag=0
                if lib.isDeleted(c):
                    flag=1
                f.write(str(c.coords[0]) + 'x' + str(c.coords[1]) + ' %d' %(flag) + '\n')
        f.close()
        #}}} 


# 'private' Utility methods
def get_vertsOfCell(cell):
    #{{{
    dimOfContainingComplex = lib.get_dimOfContainingComplex(cell) # ACK!!
    ndimVector = np.array([0 for i in range(dimOfContainingComplex)])
    coords = cell.coords
    verts = []
    
    # here's a little grease to keep things running smoothly
    while len(coords) < 3:
        coords.append([0.0,0.0])
   
    if cell.dim == 0:
        verts = [[coords[0][0],coords[1][0],coords[2][0]]]
    if cell.dim == 1:
        verts = [[coords[0][0],coords[1][0],coords[2][0]],
                 [coords[0][1],coords[1][1],coords[2][1]]]
    if (cell.dim == 2 or cell.dim == 3):
        v0 = (coords[0][0],coords[1][0],coords[2][0])
        v1 = (coords[0][0],coords[1][1],coords[2][0])
        v2 = (coords[0][0],coords[1][1],coords[2][1])
        v3 = (coords[0][0],coords[1][0],coords[2][1])
        v4 = (coords[0][1],coords[1][0],coords[2][0])
        v5 = (coords[0][1],coords[1][1],coords[2][0])
        v6 = (coords[0][1],coords[1][1],coords[2][1])
        v7 = (coords[0][1],coords[1][0],coords[2][1])
        v_list = [v0,v1,v2,v3,v4,v5,v6,v7]     
        orderedDict = OrderedDict()
        for v in v_list:
            orderedDict[v] = True

        for v in orderedDict:
            verts.append(v)
         
    return verts
    #}}}

def plot_legend(ax=None, legend_u_val=False,legend_v_val=False,
                mixed_bd_color=None,loc='center',mode='expand',bbox=None,
                legend_list = ['$u < 0$; ', '$u > 0$; ', '$u$ monotone']):
        
    if mixed_bd_color:
        mixed_boundary_color = mixed_bd_color

    if legend_u_val:
        print('plotting legend on u')
        ## some legend magic:
        u_neg_proxy     = plt.Rectangle((0, 0), 1, 1, fc=del_color,alpha=0.5)
        u_pos_proxy     = plt.Rectangle((0, 0), 1, 1, fc=undel_color,alpha=0.5)
        mixed_bd_proxy  = plt.Rectangle((0, 0), 1, 1, fc=mixed_boundary_color,alpha=0.5)
        ax.legend([u_neg_proxy, u_pos_proxy, mixed_bd_proxy],
                  ['$u < 0$; ', '$u > 0$; ', '$u$ monotone'],
                  bbox_to_anchor=bbox,mode=mode,loc=loc,
                  ncol=3,borderaxespad=0.0,fontsize=25)                     
        
    if legend_v_val:
        print('plotting legend on v')
        ## some legend magic:
        u_neg_proxy  = plt.Rectangle((0, 0), 1, 1, fc=del_color,alpha=0.5)
        u_pos_proxy  = plt.Rectangle((0, 0), 1, 1, fc=undel_color,alpha=0.5)
        v_pos_proxy  = plt.Rectangle((0, 0), 1, 1, fc='g',alpha=0.5)
        ax.legend([u_neg_proxy,u_pos_proxy, v_pos_proxy],
                  ['$u < 0$; ','$u > 0$; ','$u$ mono. & $v > 0$'],
                  bbox_to_anchor=bbox,mode='expand',
                  loc=loc,ncol=3,borderaxespad=0.0,fontsize=25)


def saveFig(fig, filename, path=os.path.expandvars('$HOME/temp/'), 
                ext='pdf', dpi=100,toStdOut=False):
    #{{{
    print("path: ", path)
    print("filename: ", filename)
    name = '%s/%s.%s' %(path,filename,ext)
    print("name: ", name)
    fig.savefig(name,dpi=dpi,transparent='True')
    if toStdOut:
        print('\t\tfrom CellComplex, cellComplex_python.py:\nsaved %s to %s' %(filename,path))
    #}}}

# ctypes overhead 
#{{{

class Two_Cell_struct (ctypes.Structure):
    _fields_ = [ ("cell_one", Cell), ("cell_two", Cell)  ]


lib.get_userDefinedFunction.restype = ctypes.c_void_p
lib.get_vector_field.restype = ctypes.c_void_p
lib.get_exit_set_function.restype = ctypes.c_void_p
lib.get_tangency_test_function.restype = ctypes.c_void_p
#lib.get_energy_density_function.restype = ctypes.c_void_p
lib.get_block_parameterization.restype = ctypes.c_void_p

lib.eval_vector_field.argtypes = [ctypes.c_double,ctypes.c_double,ctypes.c_double]
lib.eval_vector_field.restype  = ctypes.POINTER(ctypes.c_double*3)

lib.eval_exit_set.argtypes = [ctypes.c_double,ctypes.c_double]
lib.eval_exit_set.restype  = ctypes.c_double

lib.eval_autograd.restype  = ctypes.POINTER(ctypes.c_double*6) 

lib.eval_tangency_test.argtypes = [ctypes.c_double,ctypes.c_double]
lib.eval_tangency_test.restype  = ctypes.c_double

#lib.eval_energy_density.argtypes = [ctypes.c_double,ctypes.c_double]
#lib.eval_energy_density.restype  = ctypes.c_double

lib.eval_block_parameterization.argtypes = [ctypes.c_double,ctypes.c_double]
lib.eval_block_parameterization.restype  = ctypes.POINTER(ctypes.c_double*3)   

lib.pointEvalUserDefinedFunction_3d.argtypes = [ctypes.c_double,ctypes.c_double,ctypes.c_double]
lib.pointEvalUserDefinedFunction_3d.restype  = ctypes.c_double

lib.pointEvalUserDefinedFunction_2d.argtypes = [ctypes.c_double,ctypes.c_double]
lib.pointEvalUserDefinedFunction_2d.restype  = ctypes.c_double

lib.validateDomainByValidatingTopDimensionalCells.restype = CellComplex
lib.computeNodalDomain.restype = CellComplex

lib.computePersistence.restype = CellComplex

lib.is_positive_on_entire_cell.restype  = ctypes.c_bool


#{{{  Cell
lib.dim_cell.restype  = ctypes.c_int
lib.dim_cell.argtypes = [Cell]

lib.get_coordsOfCellAtLoc.restype  = ctypes.POINTER(ctypes.c_double*2)
lib.get_coordsOfCellAtLoc.argtypes = [Cell, ctypes.c_int]

#lib.computeWholeBoundaryOfCube.restype  = Cell
#lib.computeWholeBoundaryOfCube.argtypes = [Cell]
#}}}

lib.dim_cellComplex.restype  = ctypes.c_int
lib.dim_cellComplex.argtypes = [CellComplex]

lib.new_cellComplex.restype  = CellComplex
lib.new_cellComplex.argtypes = [ctypes.POINTER(ctypes.c_double * 2), ctypes.c_size_t]

lib.new_cell.restype  = Cell
# return type defined in Cell __init__

lib.size_wholeBoundaryOfCube.restype = ctypes.c_int
lib.size_wholeBoundaryOfCube.argtypes = [Cell]

lib.size_elementsAtSpecifiedDim.restype  = ctypes.c_int
lib.size_elementsAtSpecifiedDim.argtypes = [CellComplex, ctypes.c_int]

#lib.get_elementsAtSpecifiedDim.restype  = Cell
#lib.get_elementsAtSpecifiedDim.argtypes = [CellComplex, ctypes.c_int]

lib.get_boundaryElementAtPosition.restype = Cell
lib.get_boundaryElementAtPosition.argtypes = [Cell,ctypes.c_int]
    
lib.get_elementAtSpecifiedDimAndLoc.restype  = Cell 
lib.get_elementAtSpecifiedDimAndLoc.argtypes = [CellComplex, ctypes.c_int, ctypes.c_int]

lib.divideTheCubeByGivenFactor.restype  = CellComplex
lib.divideTheCubeByGivenFactor.argtypes = [CellComplex, Cell, ctypes.c_int]
 
lib.divideCubeAlongLongestAxis.restype  = Two_Cell_struct
lib.divideCubeAlongLongestAxis.argtypes = [CellComplex, Cell]

lib.divideCube.restype  = Two_Cell_struct
lib.divideCube.argtypes = [CellComplex, Cell, ctypes.c_int]
        
lib.get_dimOfContainingComplex.restype  = ctypes.c_int
lib.get_dimOfContainingComplex.argtypes = [Cell]

lib.delet.restype  = Cell
lib.delet.argtypes = [Cell]

lib.undelet.restype  = Cell
lib.undelet.argtypes = [Cell]

lib.isDeleted.restype  = ctypes.c_bool
lib.isDeleted.argtypes = [Cell]
#}}}
