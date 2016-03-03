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

#nest.SetStatus(neuron, {"I_e": 376.0})


noise = nest.Create("poisson_generator", 2)
nest.SetStatus(noise, [{"rate": 80000.0}, {"rate": 15000.0}])

syn_dict_ex = {"weight": 1.2}
syn_dict_in = {"weight": -2.0}
nest.Connect([noise[0]], neuron, syn_spec=syn_dict_ex)
nest.Connect([noise[1]], neuron, syn_spec=syn_dict_in)

multimeter = nest.Create("multimeter")
nest.SetStatus(multimeter, {"withtime":True, "record_from":["V_m"]})

#nest.SetStatus(neuron, {"V_m": 376.0})
#print (nest.GetStatus(neuron, "V_m"))

spikedetector = nest.Create("spike_detector",
                params={"withgid": True, "withtime": True})

nest.Connect(multimeter, neuron)
nest.Connect(neuron, spikedetector)

nest.Simulate(1000.0)

dmm = nest.GetStatus(multimeter)[0]
Vms = dmm["events"]["V_m"]
ts = dmm["events"]["times"]

pylab.figure(1)
pylab.plot(ts, Vms)

dSD = nest.GetStatus(spikedetector,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
pylab.figure(2)
pylab.plot(ts, evs, ".")
pylab.show()