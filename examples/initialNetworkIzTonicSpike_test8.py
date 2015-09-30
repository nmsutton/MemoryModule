import pylab
import nest
import math as math

'''
Create objects to run experiment with

Method of connecting populations:
http://www.nest-simulator.org/introduction-to-pynest/part-2-populations-of-neurons/
'''
multimeter = nest.Create("multimeter",10)
nest.SetStatus(multimeter, {"withtime":True, "record_from":["V_m"]})
multimeter2 = nest.Create("multimeter")
nest.SetStatus(multimeter2, {"withtime":True, "record_from":["V_m"]})
spikedetector = nest.Create("spike_detector",
                params={"withgid": True, "withtime": True})
spikedetector2 = nest.Create("spike_detector",
                params={"withgid": True, "withtime": True})

'''noise = nest.Create("poisson_generator", 2)
nest.SetStatus(noise, [{"rate": 80000.0}, {"rate": 15000.0}])'''

# values for neurons taken from http://neuralensemble.org/docs/PyNN/examples/Izhikevich.html?highlight=izhikevich
pop1 = nest.Create("izhikevich",100,{'V_m':-70.0,'I_e':-160.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})
pop2 = nest.Create("izhikevich",100,{'V_m':-70.0,'I_e':-180.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})

'''
Form connections between objects and run sim
'''

syn_weight = 100.0

def createSyn(input_layer, output_layer, fire_rate_ratio, syn_weight):
	'''
		Note: later uneven numbers of neurons in layers
		could be added but for now using even.

		Ratio of 1.0 creates 50% ex and 50% inh
		2.0 creates 66% ex and 33% inh
		0.5 creates 33% ex and 66% inh

		TODO: check if ratio calc works exactly right

		"fixed_total_number. Here n connections are created ... from the 
		populations pre and ... post."

		TODO: for now synapses are one-to-one to control ratio of responses.
		In the future more e.g. one-to-many should be made while controlling
		activity between layers

		Note: It is needed that the number of exhitatory connections totally
		control the amount of firing fore each next layer, no origional firing
		occurs in any layer but the first.
	'''
	len_in_layer = len(input_layer)
	len_out_layer = len(output_layer)
	times_greater_ratio = math.ceil(fire_rate_ratio)
	syn_dict = {"weight": syn_weight}

	for time_greater in range(times_greater_ratio):
		adjusted_delay = 0.1 + (60.0 * time_greater)
		adjusted_conn_total = len_in_layer
		if (time_greater==(times_greater_ratio-1)): 
			adjusted_conn_total = math.floor(len_in_layer*(fire_rate_ratio-(times_greater_ratio-1)))

		syn_dict = {"weight": syn_weight, "delay":adjusted_delay} 
		for i in range(adjusted_conn_total):
			nest.Connect([input_layer[i]], [output_layer[i]], "one_to_one", syn_dict)#syn_spec=syn_weight)

createSyn(pop1,pop2,0.928, syn_weight)
#createSyn(pop1,pop2,1.000, syn_weight)
#createSyn(pop1,pop2,2.164, syn_weight)
#createSyn(pop1,pop2,3.0, syn_weight)

nest.Connect(multimeter, pop1)
nest.Connect(multimeter2, pop2)
nest.Connect(pop1, spikedetector)
nest.Connect(pop2, spikedetector2)

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
pylab.figure(1)
pylab.plot(ts, Vms)

pylab.figure(2)
pylab.plot(ts2, Vms2)

dSD = nest.GetStatus(spikedetector,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
print ('number of spikes')
print (len(ts))
print(sum(ts>800))
'''pylab.figure(3)
pylab.plot(ts, evs, ".")'''

dSD_b = nest.GetStatus(spikedetector2,keys='events')[0]
evs_b = dSD_b["senders"]
ts_b = dSD_b["times"]
print ('number of spikes')
print (len(ts_b))
print(sum(ts_b>800))
#pylab.figure(4)
#pylab.plot(ts_b, evs_b, ".")

pylab.show()