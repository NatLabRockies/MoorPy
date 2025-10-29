
import pytest

from numpy.testing import assert_allclose

import numpy as np
import moorpy as mp
import matplotlib.pyplot as plt
from moorpy.helpers import lines2ss
import os

def test_seabed_load():
    # check bathymetry grid loads in correctly
    direct = os.path.dirname(os.path.realpath(__file__))
    ms = mp.System(file=os.path.join(direct,'MoordynSemiTautUpd.dat'), 
                   bathymetry=os.path.join(direct,'bathymetry200m_sample.txt'))
    ms.initialize()
    ms.solveEquilibrium()
    assert len(ms.bathGrid_Xs)==2 
    assert_allclose(ms.bathGrid,[[200, 210.],[210., 220.]])
    ms.plot()
    # convert to subsystems and check again
    ms2 = lines2ss(ms)
    ms2.initialize()
    ms2.solveEquilibrium()
    assert len(ms2.bathGrid_Xs)==2
    assert_allclose(ms.bathGrid,[[200, 210.],[210., 220.]])
    assert len(ms2.lineList[0].bathGrid_Xs)==2
    assert_allclose(ms2.lineList[0].bathGrid,[[218.13014286, 211.59078571]])
    

def test_seabed_lay():
    # check mooring line is laying on the seabed considering bathymetry
    direct = os.path.dirname(os.path.realpath(__file__))
    ms = mp.System(file=os.path.join(direct,'MoordynSemiTautUpd.dat'), 
                   bathymetry=os.path.join(direct,'bathymetry200m_sample.txt'))
    ms.initialize()
    ms.solveEquilibrium()
    # check mooring line Z depth is properly evolving
    assert ms.lineList[0].Zs[0]<ms.lineList[0].Zs[1] 
    
    # convert to subsystems and check again
    ms2 = lines2ss(ms)
    ms2.initialize()
    ms2.solveEquilibrium()
    # check mooring line Z depth is properly evolving
    assert ms.lineList[0].lineList[0].Zs[0]==pytest.approx(-218.13014286)
    assert ms.lineList[0].lineList[0].Zs[0]<ms.lineList[0].lineList[0].Zs[1] 


    