'''
Weakly adapting PYR model
@author: Ferguson et al. (2014) F1000Res. 
'''
from brian import *

defaultclock.dt = 0.02*ms

#Weakly adapting PYR parameters for model 1  
C=300 * pF
vr=-61.8 * mV
vpeak=22.6 * mV
c=-65.8 * mV
klow=0.5 * nS/mV
khigh=3.3  * nS/mV
a= 0.001 /ms   #0.00008 /ms
d=5 * pA
vt=-57.0 *mV
b=3 * nS

N=1   #number of cells
mean_Iapp=100 #mean Iapplied input (pA) 
Ishift_raw=-45  #Ishift (pA)

time=0

#cell eqns
pyr_eqs = """
Iext  : amp
Ishift : amp
k=(v<vt)*klow+(v>=vt)*khigh : (siemens/volt)
du/dt = a*(b*(v-vr)-u)            : amp
dv/dt = (k*(v-vr)*(v-vt)+Ishift+Iext -u)/C : volt
"""

#define neuron group
PYR = NeuronGroup(N, model=pyr_eqs, reset ="v = c; u += d" , threshold="v>=vpeak")

#set excitatory drive 
PYR.Iext = mean_Iapp*pA

#set Ishift 
PYR.Ishift = Ishift_raw*pA

#set initial conditions for each neuron
PYR.v = rand(len(PYR))*0.01 -0.065

#record all spike times for the neuron group
PYR_v = StateMonitor(PYR, 'v', record=True)

#run for x seconds of simulated time
duration = 1 * second  # 0.01 * second

#create arrays to record network output
times = zeros(duration/defaultclock.dt +1)
voltage = zeros(duration/defaultclock.dt +1)

net =Network(PYR,PYR_v) 
net.run(duration)


####make voltage plot####
plot(PYR_v.times,PYR_v[0]/mV)
xlabel("Time (s)")
ylabel("Membrane Potential (mV)")
title('Weakly adapting PYR model 1 with %d pA input'%(mean_Iapp))
#title('Weakly adapting PYR model 2 with %d pA input'%(mean_Iapp))  
show()
