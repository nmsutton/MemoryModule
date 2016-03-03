# Tutorial2.py
# Author : Jose Guzman
# Mail : jose.guzman __at__ physiologie.uni-freiburg.de

#! /usr/bin/env python

import neuron
import nrn
from neuron import h

# create soma object
soma = nrn.Section()

# soma attributes
soma.L = 30
soma.nseg = 3
soma.diam = 30

h('''
	NDEND = 2
	create soma, basal[NDEND], axon
	access basal[1]
	''')

sec = h.cas()
print (sec(), sec.name())