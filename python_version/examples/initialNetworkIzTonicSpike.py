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

# values for neurons taken from http://neuralensemble.org/docs/PyNN/examples/Izhikevich.html?highlight=izhikevich
pop1 = nest.Create("izhikevich",1,{'V_m':-70.0,'I_e':14.0,'a':0.02,'b':0.2,'c':-65.0,'d':6.0})
pop2 = nest.Create("izhikevich",1,{'V_m':-70.0,'I_e':14.0,'a':0.02,'b':0.2,'c':-65.0,'d':6.0})
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
nest.SetStatus(pop1, {"I_e": 10.0})

nest.Connect(pop1, pop2, syn_spec = {"weight":20.0})

nest.Connect(multimeter, pop1)
nest.Connect(multimeter2, pop2)
nest.Connect(pop1, spikedetector)

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

'''pylab.figure(2)
pylab.plot(ts2, Vms2)

dSD = nest.GetStatus(spikedetector,keys='events')[0]
evs = dSD["senders"]
ts = dSD["times"]
pylab.figure(3)
pylab.plot(ts, evs, ".")'''
pylab.show()