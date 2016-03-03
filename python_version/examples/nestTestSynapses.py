#import sys, os
#sys.path.insert(0, '/media/nmsutton/Ext3Drive/General/NEST/NEST/lib64/python3.4/site-packages')

#print (os.path.dirname(sys.executable))
#print (os.environ)

import pylab
import nest

print ("test")

neuron = nest.Create("iaf_neuron")
neuron2 = nest.Create("iaf_neuron")
neuron3 = nest.Create("iaf_neuron")

multimeter = nest.Create("multimeter")
nest.SetStatus(multimeter, {"withtime":True, "record_from":["V_m"]})
multimeter2 = nest.Create("multimeter")
nest.SetStatus(multimeter2, {"withtime":True, "record_from":["V_m"]})

nest.SetStatus(neuron, {"I_e": 376.0})
nest.SetStatus(neuron2, {"I_e": 326.0})
#nest.SetStatus(neuron2, {"V_m": 376.0})

nest.Connect(neuron, neuron2, syn_spec = {"weight":-150.0})

#nest.SetStatus(neuron, {"V_m": 376.0})
#print (nest.GetStatus(neuron, "V_m"))

nest.Connect(multimeter, neuron)
nest.Connect(multimeter2, neuron2)

spikedetector = nest.Create("spike_detector",
                params={"withgid": True, "withtime": True})
nest.Connect(neuron, spikedetector)

nest.Simulate(1000.0)

dmm = nest.GetStatus(multimeter)[0]
Vms = dmm["events"]["V_m"]
ts = dmm["events"]["times"]

pylab.figure(1)
pylab.plot(ts, Vms)

dmb = nest.GetStatus(multimeter2)[0]
Vmsb = dmb["events"]["V_m"]
tsb = dmb["events"]["times"]

pylab.figure(2)
pylab.plot(tsb, Vmsb)

'''dSD = nest.GetStatus(spikedetector,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
pylab.figure(2)
pylab.plot(ts, evs, ".")'''

pylab.show()