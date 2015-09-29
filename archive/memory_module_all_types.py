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
spikedetector = nest.Create("spike_detector",
                params={"withgid": True, "withtime": True})

'''noise = nest.Create("poisson_generator", 2)
nest.SetStatus(noise, [{"rate": 80000.0}, {"rate": 15000.0}])'''

e_c_2_layer = nest.Create("izhikevich",50,{'V_m':-70.0,'I_e':-15.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})
e_c_3_layer = nest.Create("izhikevich",90,{'V_m':-70.0,'I_e':-15.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})
e_c_5_layer = nest.Create("izhikevich",80,{'V_m':-70.0,'I_e':-15.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})
c_a_1_layer = nest.Create("izhikevich",340,{'V_m':-70.0,'I_e':-15.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})
c_a_3_layer = nest.Create("izhikevich",100,{'V_m':-70.0,'I_e':-15.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})
d_g_layer = nest.Create("izhikevich",12,{'V_m':-70.0,'I_e':-15.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})

'''
Form connections between neurons and run sim

NOTE: I may need to split the neurons into Ex and In 
groups in layers for connections

With a number of neuron mismatch between layers
how is that processed in connections?
'''
syn_dict_ex = {"weight": 1.2}
syn_dict_in = {"weight": -2.0}
nest.Connect(e_c_2_layer, e_c_3_layer, "all_to_all", syn_spec=syn_dict_ex)
nest.Connect(e_c_3_layer, e_c_5_layer, "all_to_all", syn_spec=syn_dict_ex)
nest.Connect(e_c_5_layer, c_a_1_layer, "all_to_all", syn_spec=syn_dict_ex)
nest.Connect(c_a_1_layer, c_a_3_layer, "all_to_all", syn_spec=syn_dict_ex)
nest.Connect(c_a_3_layer, d_g_layer, "all_to_all", syn_spec=syn_dict_ex)

nest.SetStatus(e_c_2_layer, {"I_e": 10.0})

nest.Connect(multimeter, e_c_2_layer)
nest.Connect(multimeter2, d_g_layer)
nest.Connect(e_c_2_layer, spikedetector)

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