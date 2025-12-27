# -*- coding: utf-8 -*-
"""
Tests nonlinear stiffness functionality (alternative to Catenary)
"""

import moorpy as mp
import numpy as np
from numpy.testing import assert_allclose

def test_nonlinear():
    # create a simple system
    ms = mp.System(depth=600)
    ms.setLineType(220,'polyester')
    lt0 = ms.lineTypes['0']
    lt = dict(m=lt0['m'],
              d_nom=lt0['d_nom'],
              d_vol=lt0['d_vol'],
              material=lt0['material'],
              )
    
    # add some strains and tensions
    Str = np.array([.001,.025,.04,.065,.1,.125]) # [-] (L_new/L-1)
    Ten = np.array([1.25e6,1.5e6,3.3e6,5.8e6,8.5e6,1.3e7]) # [N]
    
    # create a line type
    ms.setLineType(lineType=lt,
                   name="1",
                   Str=Str,
                   Ten=Ten,
                   )
    
    # set up a single line attached to a single body
    len_0 = 1190
    ms.addPoint(1,[-1131.37,0,-600])
    ms.addPoint(1,[-58,0,-14])
    ms.addBody(-1,[0,0,0,0,0,0])
    ms.bodyList[0].attachPoint(2,[-58,0,-14])
    ms.addLine(len_0,'1',pointA=1,pointB=2)
    
    # solve equilibrium
    ms.initialize()
    ms.solveEquilibrium()
    loi = ms.lineList[0] # line of interest
    len_init = np.linalg.norm(loi.rB-loi.rA) # initial length
    
    # calculate the tension at end B based on the str-ten relationship & line props
    str_init = (len_init/len_0)-1
    ten_init = np.interp(str_init,Str,Ten)
    ten_cosine = np.linalg.norm(loi.rB[:2]-loi.rA[:2])/len_init
    ten_sine = (loi.rB[2]-loi.rA[2])/len_init
    ten_init_B = np.linalg.norm([-ten_init*ten_cosine,-ten_init*ten_sine-.5*loi.type['w']*loi.L])
    
    # assert the calculated tension matches the TB
    assert_allclose(ten_init_B,loi.TB,rtol=.001)
    
    # move the body
    ms.bodyList[0].r6=[25,0,0,0,0,0]
    # re-calc equilibrium
    ms.initialize()
    ms.solveEquilibrium()
    
    # calculate the new tension at end B
    len_new = np.linalg.norm(loi.rB-loi.rA)
    str_new = (len_new/len_0)-1
    ten_new=np.interp(str_new,Str,Ten)
    ten_cosine = np.linalg.norm(loi.rB[:2]-loi.rA[:2])/len_new
    ten_sine = (loi.rB[2]-loi.rA[2])/len_new
    ten_new_B = np.linalg.norm([-ten_new*ten_cosine,-ten_new*ten_sine-.5*loi.type['w']*loi.L])
    
    # assert the calculated tension matchse the TB
    assert_allclose(ten_new_B,loi.TB,rtol=.001)