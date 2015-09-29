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

'''noise = nest.Create("poisson_generator", 2)
nest.SetStatus(noise, [{"rate": 80000.0}, {"rate": 15000.0}])'''

# values for neurons taken from http://neuralensemble.org/docs/PyNN/examples/Izhikevich.html?highlight=izhikevich
pop1 = nest.Create("izhikevich",10,{'V_m':-70.0,'I_e':18.0,'a':0.005,'b':0.2,'c':-65.0,'d':6.0})
pop2 = nest.Create("izhikevich",10,{'V_m':-70.0,'I_e':4.0,'a':0.02,'b':0.2,'c':-65.0,'d':6.0})
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

def createSyn(input_layer, output_layer, fire_rate_ratio):
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
	syn_dict = {"weight": syn_weight} 

	perc_conn = fire_rate_ratio / (fire_rate_ratio+1)

	conn_total = math.floor(len_in_layer*perc_conn)
	#print (conn_total)
	
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


createSyn(pop1,pop2,0.928)

nest.Connect(multimeter, pop1)
nest.Connect(multimeter2, pop2)
nest.Connect(pop2, spikedetector)

nest.Simulate(350.0)

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
pylab.figure(3)
pylab.plot(ts, evs, ".")
pylab.show()