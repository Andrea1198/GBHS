#!/usr/bin/env python
__author__ = "Tommaso Chiarotti"
__license__ = "see LICENSE.txt file."
__version__ = "0"
__email__ = "tommaso.chiarotti@epfl.ch"
"""
This program import module_graph.py to start plotting procedure
"""

import module_graph as mgraph

dirfile   = 'RPA_test.save'
listiter  = [1] 
listdelta = [0.01]
listpps   = [0.05]
vstart    = 1000
vend      = 1100
vstep     = 100
xstart    = 1000 
xend      = 1100
xstep     = 100
xylabel   = [[None,None],[None,None],[None,None],[None,None]]
xylim     = [[None,None],[None,None],[None,None],[None,None]]

plotter = mgraph.plot(dirfile,listiter,listdelta,listpps,vstart,vend,vstep,xstart,xend,xstep)

#plotter.sigma_and_greenf(xylabel,xylim)

#plotter.sc_sigma_and_greenf(xylabel,xylim)

#plotter.occ_number(xylabel,xylim)

#plotter.sc_occ_number(xylabel,xylim)

#plotter.pol_and_wscreen(xylabel,xylim)

#plotter.sc_pol_and_wscreen(xylabel,xylim)
