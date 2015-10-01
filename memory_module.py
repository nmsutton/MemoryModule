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
import math as math
import numpy as np

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

e_c_3_layer = nest.Create("izhikevich",500,{'V_m':-70.0,'I_e':-160.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})
e_c_5_layer = nest.Create("izhikevich",500,{'V_m':-70.0,'I_e':-180.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})
c_a_1_layer = nest.Create("izhikevich",500,{'V_m':-70.0,'I_e':-180.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})

'''
Form connections between neurons and run sim

NOTE: I may need to split the neurons into Ex and In 
groups in layers for connections

With a number of neuron mismatch between layers
how is that processed in connections?
'''

'''
	Synapses
'''
syn_weight = 300.0

def createSyn(input_layer, output_layer, fire_rate_ratio, syn_weight, neuron_range):
	'''
		neuron_range: min_index, max_index, total_neurons

		Note: later uneven numbers of neurons in layers
		could be added but for now using even.

		Ratio of 1.0 creates 50% ex and 50% inh
		2.0 creates 66% ex and 33% inh
		0.5 creates 33% ex and 66% inh

		TODO: check if ratio calc works exactly right

		TODO: for now synapses are one-to-one to control ratio of responses.
		In the future more e.g. one-to-many should be made while controlling
		activity between layers

		Note: It is needed that the number of exhitatory connections totally
		control the amount of firing fore each next layer, no origional firing
		occurs in any layer but the first.
	'''

	times_greater_ratio = math.ceil(fire_rate_ratio)
	syn_dict = {"weight": syn_weight}
	min_index = neuron_range[0]
	max_index = neuron_range[1]
	total_neurons = neuron_range[2]
	total_range = range(math.floor((max_index-min_index)*total_neurons))
	total_range = np.array(total_range) + math.floor(min_index*total_neurons)
	len_in_layer = len(total_range)
	#len_out_layer = len(output_layer)

	for time_greater in range(times_greater_ratio):
		adjusted_delay = 0.1 + (20.0 * time_greater)
		adjusted_conn_total = len_in_layer
		if (time_greater==(times_greater_ratio-1)): 
			adjusted_conn_total = math.floor(len_in_layer*(fire_rate_ratio-(times_greater_ratio-1)))

		syn_dict = {"weight": syn_weight, "delay":adjusted_delay} 
		for i in range(adjusted_conn_total):
			#print(total_range[0])
			#print(i)
			n_i = total_range[i]
			nest.Connect([input_layer[n_i]], [output_layer[n_i]], "one_to_one", syn_dict)

createSyn(e_c_3_layer,e_c_5_layer,1.4917, syn_weight, [0.0, 0.233, len(e_c_3_layer)])
createSyn(e_c_3_layer,e_c_5_layer,2.2081, syn_weight, [0.233, 0.567, len(e_c_3_layer)])
createSyn(e_c_3_layer,e_c_5_layer,0.6152, syn_weight, [0.567, 0.667, len(e_c_3_layer)])
createSyn(e_c_3_layer,e_c_5_layer,0.3024, syn_weight, [0.667, 1.0, len(e_c_3_layer)])

createSyn(e_c_5_layer,c_a_1_layer,6.8897, syn_weight, [0.0, 0.233, len(e_c_3_layer)])
createSyn(e_c_5_layer,c_a_1_layer,4.6546, syn_weight, [0.233, 0.367, len(e_c_3_layer)])
createSyn(e_c_5_layer,c_a_1_layer,1.6016, syn_weight, [0.367, 0.567, len(e_c_3_layer)])
createSyn(e_c_5_layer,c_a_1_layer,5.7480, syn_weight, [0.567, 0.733, len(e_c_3_layer)])
createSyn(e_c_5_layer,c_a_1_layer,7.6722, syn_weight, [0.733, 1.0, len(e_c_3_layer)])

#createSyn(e_c_3_layer,e_c_5_layer,0.928, syn_weight)
#createSyn(e_c_3_layer,e_c_5_layer,0.3, syn_weight, [1, 7, len(e_c_3_layer)])
#createSyn(e_c_3_layer,e_c_5_layer,0.3, syn_weight, [0.0, 0.5, len(e_c_3_layer)])
#createSyn(e_c_3_layer,e_c_5_layer,0.3, syn_weight, [0.50, 1.0, len(e_c_3_layer)])
# 2.164/0.928=2.332
#createSyn(e_c_5_layer,c_a_1_layer,2.332, syn_weight)
#createSyn(e_c_5_layer,c_a_1_layer,8.0, syn_weight, [1, 7, len(e_c_3_layer)])
#createSyn(e_c_5_layer,c_a_1_layer,8.0, syn_weight, [0.0, 1.0, len(e_c_3_layer)])

nest.Connect(multimeter, e_c_3_layer)
nest.Connect(multimeter2, c_a_1_layer)

nest.Connect(e_c_3_layer, spikedetector_e_c_3)
nest.Connect(e_c_5_layer, spikedetector_e_c_5)
nest.Connect(c_a_1_layer, spikedetector_c_a_1)

'''
	NOTE: filtering of spike counts after a certain
	time happens later and therefore only a portion
	of sim time is counted.
'''
nest.Simulate(2000.0)

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
print(sum(ts>800))
#pylab.figure(2)
#pylab.plot(ts, evs, ".")

dSD = nest.GetStatus(spikedetector_e_c_5,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
print ('number of spikes')
print(sum(ts>800))
#pylab.figure(3)
#pylab.plot(ts, evs, ".")

dSD = nest.GetStatus(spikedetector_c_a_1,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
print ('number of spikes')
print(sum(ts>800))
#pylab.figure(4)
#pylab.plot(ts, evs, ".")

pylab.show()