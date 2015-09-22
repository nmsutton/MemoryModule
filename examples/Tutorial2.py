# Tutorial2.py
# Author : Jose Guzman
# Mail : jose.guzman __at__ physiologie.uni-freiburg.de

#! /usr/bin/env python

import neuron
import nrn
from neuron import h

# load Neuron graphical user interface
#if not ( h('load_file("nrngui.hoc")')):
#    print "Error, cannot open NEURON gui"

# create soma object
soma = nrn.Section()

# soma attributes
soma.L = 30
soma.nseg = 3
soma.diam = 30

soma.insert('hh')

#stimulation

stimulator = h.IClamp(0.5,soma)

stimulator.delay = 0.1
stimulator.dur = 0.8
stimulator.amp = 1.2

# test Vmb without perturbation
#print "The initial membrane potential is %6.4f" %(soma(.5).v)
print (soma(.5).v)


# Run simulation
print (soma.v)

# set tstop at 1.5
tstop = 1.5

neuron.run(tstop)

print (soma.v)

tstop = 5

neuron.run(tstop)

print (soma.v)

def get_vmb(time):
    """ Returns the voltage at the soma at a given time"""

    neuron.run(time)

    if h('run()'):
        return soma(.5).v
    else:
        return 0
