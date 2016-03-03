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
pop1 = nest.Create("izhikevich",100,{'V_m':-70.0,'I_e':14.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})
pop2 = nest.Create("izhikevich",100,{'V_m':-70.0,'I_e':14.0,'a':0.0012,'b':3.0,'c':-68.5,'d':10.0})
#pop1 = nest.Create("izhikevich",{'a':0.02,'b':0.2,'d':6.0})
#pop2 = nest.Create("izhikevich",{'a':0.02,'b':0.2,'d':6.0})
#pop1 = nest.Create("izhikevich")
#pop2 = nest.Create("izhikevich")

'''
Form connections between objects and run sim
'''
'''syn_dict_ex = {"weight": 1.2}
syn_dict_in = {"weight": -2.0}
nest.Connect([noise[0]], pop1, syn_spec=syn_dict_ex)
nest.Connect([noise[1]], pop1, syn_spec=syn_dict_in)'''
#nest.SetStatus(pop1, {"I_e": 376.0})
#nest.SetStatus(pop1, {"I_e": 10.0})

#nest.Connect(pop1, pop2, syn_spec = {"weight":-10.0})

# find number of neurons in layer
print('len pop1:')
print(len(pop1))

syn_weight = 50.0

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
	max_syn_ratio = 2.5
	syn_multipler = 2.0

	perc_conn = fire_rate_ratio

	#if fire_rate_ratio >= 1.0:
		#perc_conn = fire_rate_ratio / (fire_rate_ratio+1)
	#elif fire_rate_ratio < 1.0:
	if fire_rate_ratio < 1.0:
		#perc_conn = (1/fire_rate_ratio) / ((1/fire_rate_ratio)+1)
		prec_conn = 1/perc_conn
		syn_weight = syn_weight * -1.0
	perc_conn = 1 - perc_conn
	print (perc_conn)

	syn_dict = {"weight": syn_weight} 

	conn_total = math.floor(len_in_layer*(syn_multipler*perc_conn/max_syn_ratio))
	print (conn_total)
	
	'''conn_dict = {"rule": "fixed_outdegree", "outdegree":conn_total, "autapses": False, "multapses": False}
	syn_dict = {"weight": syn_weight} 
	nest.Connect(pop1, pop2, conn_dict, syn_dict)'''

	'''for i in range(numb_ex):
		nest.Connect(input_layer[i], output_layer[i], syn_spec=syn_dict_ex)
	for i in range(numb_inh):
		nest.Connect(input_layer[i], output_layer[i], syn_spec=syn_dict_in)'''

	'''print (input_layer[0])
	print (output_layer[0])
	print (input_layer)
	print (output_layer)'''

	for i in range(conn_total):
		print(i)
		nest.Connect([input_layer[i]], [output_layer[i]], "one_to_one", syn_dict)#syn_spec=syn_weight)
		print ('node')
		print (i)


createSyn(pop1,pop2,0.928, syn_weight)

nest.Connect(multimeter, pop1)
nest.Connect(multimeter2, pop2)
nest.Connect(pop1, spikedetector)
nest.Connect(pop2, spikedetector2)

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
'''pylab.figure(1)
pylab.plot(ts, Vms)

pylab.figure(2)
pylab.plot(ts2, Vms2)'''

dSD = nest.GetStatus(spikedetector,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
print ('number of spikes')
print (len(ts))
'''pylab.figure(3)
pylab.plot(ts, evs, ".")'''

dSD_b = nest.GetStatus(spikedetector2,keys='events')[0]
evs_b = dSD_b["senders"]
ts_b = dSD_b["times"]
print ('number of spikes')
print (len(ts_b))
#pylab.figure(4)
#pylab.plot(ts_b, evs_b, ".")

pylab.show()