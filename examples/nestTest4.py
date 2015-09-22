'''
Examples used from: http://www.nest-simulator.org/introduction-to-pynest/part-1-neurons-and-simple-neural-networks/
'''

import pylab
import nest

'''neuron1 = nest.Create("iaf_neuron")
#nest.SetStatus(neuron1, {"I_e": 376.0})
neuron2 = nest.Create("iaf_neuron")
multimeter = nest.Create("multimeter")
nest.SetStatus(multimeter, {"withtime":True, "record_from":["V_m"]})
multimeter2 = nest.Create("multimeter")
nest.SetStatus(multimeter2, {"withtime":True, "record_from":["V_m"]})'''
pop1 = nest.Create("iaf_neuron",10)
nest.SetStatus(pop1, {"I_e": 376.0})
pop2 = nest.Create("iaf_neuron",10)
multimeter = nest.Create("multimeter",10)
nest.SetStatus(multimeter, {"withtime":True, "record_from":["V_m"]})
multimeter2 = nest.Create("multimeter")
nest.SetStatus(multimeter2, {"withtime":True, "record_from":["V_m"]})

noise = nest.Create("poisson_generator", 2)
nest.SetStatus(noise, [{"rate": 80000.0}, {"rate": 15000.0}])

syn_dict_ex = {"weight": 1.2}
syn_dict_in = {"weight": -2.0}
nest.Connect([noise[0]], pop1, syn_spec=syn_dict_ex)
nest.Connect([noise[1]], pop1, syn_spec=syn_dict_in)

#nest.Connect(pop1, pop2, syn_spec = {"weight":20.0})
nest.Connect(pop1, pop2, "all_to_all", {"weight":20.0})
'''
NOTE: figure 2's height of voltage traces is directly effected by the weight specified above
'''
nest.Connect(multimeter, pop1)
nest.Connect(multimeter2, pop2)

spikedetector = nest.Create("spike_detector",
                params={"withgid": True, "withtime": True})

#nest.Connect(multimeter, neuron)
nest.Connect(pop1, spikedetector)

nest.Simulate(1000.0)

dmm = nest.GetStatus(multimeter)[0]
Vms = dmm["events"]["V_m"]
ts = dmm["events"]["times"]

dmm2 = nest.GetStatus(multimeter2)[0]
Vms2 = dmm2["events"]["V_m"]
ts2 = dmm2["events"]["times"]

pylab.figure(1)
pylab.plot(ts, Vms)

pylab.figure(2)
pylab.plot(ts2, Vms2)

dSD = nest.GetStatus(spikedetector,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
pylab.figure(3)
pylab.plot(ts, evs, ".")
pylab.show()