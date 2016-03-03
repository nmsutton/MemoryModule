#import sys, os
#sys.path.insert(0, '/media/nmsutton/Ext3Drive/General/NEST/NEST/lib64/python3.4/site-packages')

#print (os.path.dirname(sys.executable))
#print (os.environ)

import pylab
import nest

print ("test")

neuron = nest.Create("iaf_neuron")

nest.GetStatus(neuron)

nest.GetStatus(neuron, "I_e")
print (nest.GetStatus(neuron, ["V_reset", "V_th"]))

nest.SetStatus(neuron, {"I_e": 376.0})

multimeter = nest.Create("multimeter")
nest.SetStatus(multimeter, {"withtime":True, "record_from":["V_m"]})

nest.SetStatus(neuron, {"V_m": 376.0})
print (nest.GetStatus(neuron, "V_m"))

spikedetector = nest.Create("spike_detector",
                params={"withgid": True, "withtime": True})