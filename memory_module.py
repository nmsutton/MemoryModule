'''
Simulation Based on Hippocampus Recordings
Copyright Nate Sutton 2015

References:
Data from CRCNS.org hc3 .
Izhikevich neuron parameters from:
http://f1000research.com/articles/3-104/v1
'''

import pylab
import nest

'''
Create objects to run experiment with
'''
multimeter = nest.Create("multimeter",10)
nest.SetStatus(multimeter, {"withtime":True, "record_from":["V_m"]})
multimeter2 = nest.Create("multimeter")
nest.SetStatus(multimeter2, {"withtime":True, "record_from":["V_m"]})
spikedetector_e_c_3 = nest.Create("spike_detector", params={"withgid": True, "withtime": True})
spikedetector_e_c_5 = nest.Create("spike_detector", params={"withgid": True, "withtime": True})
spikedetector_c_a_1 = nest.Create("spike_detector", params={"withgid": True, "withtime": True})

'''noise = nest.Create("poisson_generator", 2)
nest.SetStatus(noise, [{"rate": 80000.0}, {"rate": 15000.0}])'''

e_c_3_layer = nest.Create("izhikevich",100,{'V_m':-70.0,'I_e':14.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})
e_c_5_layer = nest.Create("izhikevich",100,{'V_m':-70.0,'I_e':14.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})
c_a_1_layer = nest.Create("izhikevich",100,{'V_m':-70.0,'I_e':-160.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})

'''
Form connections between neurons and run sim

NOTE: I may need to split the neurons into Ex and In 
groups in layers for connections

With a number of neuron mismatch between layers
how is that processed in connections?
'''
# avoid default I_e values by setting them here
#nest.SetStatus(e_c_3_layer, {"I_e": -15.1})
#nest.SetStatus(e_c_5_layer, {"I_e": -50.1})
#nest.SetStatus(c_a_1_layer, {"I_e": -50.1})

syn_dict_ex = {"weight": 1.2}
syn_dict_in = {"weight": -2.0}
#syn_ec3_to_ec5 = {"weight": 0.928}
syn_ec3_to_ec5 = {"weight": -20.0}
syn_ec5_to_ca1 = {"weight": 2.164}
nest.Connect(e_c_3_layer, e_c_5_layer, "one_to_one", syn_spec=syn_ec3_to_ec5)
nest.Connect(e_c_5_layer, c_a_1_layer, "one_to_one", syn_spec=syn_ec5_to_ca1)

nest.Connect(multimeter, e_c_3_layer)
nest.Connect(multimeter2, c_a_1_layer)
nest.Connect(e_c_3_layer, spikedetector_e_c_3)
nest.Connect(e_c_5_layer, spikedetector_e_c_5)
nest.Connect(c_a_1_layer, spikedetector_c_a_1)

#nest.Simulate(350.0)
nest.Simulate(1000.0)

'''
Record activity
'''
dmm = nest.GetStatus(multimeter)[0]
Vms = dmm["events"]["V_m"]
ts = dmm["events"]["times"]

dmm2 = nest.GetStatus(multimeter2)[0]
Vms2 = dmm2["events"]["V_m"]
ts2 = dmm2["events"]["times"]

'''
Plot results
'''
#pylab.figure(1)
#pylab.plot(ts, Vms)

#pylab.figure(2)
#pylab.plot(ts2, Vms2)

dSD = nest.GetStatus(spikedetector_e_c_3,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
print ('number of spikes')
print (len(ts))
pylab.figure(2)
pylab.plot(ts, evs, ".")

dSD = nest.GetStatus(spikedetector_e_c_5,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
print ('number of spikes')
print (len(ts))
pylab.figure(3)
pylab.plot(ts, evs, ".")

dSD = nest.GetStatus(spikedetector_c_a_1,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
print ('number of spikes')
print (len(ts))
pylab.figure(4)
pylab.plot(ts, evs, ".")

pylab.show()