/*
 * Copyright (c) 2013 Regents of the University of California. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. The names of its contributors may not be used to endorse or promote
 *    products derived from this software without specific prior written
 *    permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * *********************************************************************************************** *
 * CARLsim
 * created by: 		(MDR) Micah Richert, (JN) Jayram M. Nageswaran
 * maintained by:	(MA) Mike Avery <averym@uci.edu>, (MB) Michael Beyeler <mbeyeler@uci.edu>,
 *					(KDC) Kristofor Carlson <kdcarlso@uci.edu>
 *
 * CARLsim available from http://socsci.uci.edu/~jkrichma/CARL/CARLsim/
 * Ver 12/3/2014
 */

#include <carlsim.h>
#include <stdio.h>

#if defined(WIN32) || defined(WIN64)
	#define _CRT_SECURE_NO_WARNINGS
	#include <periodic_spikegen.h>
	#include <simple_weight_tuner.h>
#endif

int main(int argc, const char* argv[]) {
	// ---------------- CONFIG STATE -------------------
	CARLsim *sim = new CARLsim("MemModGPU", GPU_MODE, USER, 0, 42);
	int neuronsPerGroup = 500;

	// input is a SpikeGenerator group that fires every 20 ms (50 Hz)
	PeriodicSpikeGenerator PSG(1.15f);
	int spike_gen=sim->createSpikeGeneratorGroup("sg", neuronsPerGroup, EXCITATORY_NEURON);
	sim->setSpikeGenerator(spike_gen, &PSG);

	//
	int e_c_3_layer1=sim->createGroup("ec3_1", ceil(neuronsPerGroup*.233), EXCITATORY_NEURON);
	int e_c_3_layer2=sim->createGroup("ec3_2", ceil(neuronsPerGroup*.133), EXCITATORY_NEURON);
	int e_c_3_layer3=sim->createGroup("ec3_3", ceil(neuronsPerGroup*.2), EXCITATORY_NEURON);
	int e_c_3_layer4=sim->createGroup("ec3_4", ceil(neuronsPerGroup*.1), EXCITATORY_NEURON);
	int e_c_3_layer5=sim->createGroup("ec3_5", ceil(neuronsPerGroup*.067), EXCITATORY_NEURON);
	int e_c_3_layer6=sim->createGroup("ec3_6", ceil(neuronsPerGroup*.267), EXCITATORY_NEURON);
	sim->setNeuronParameters(e_c_3_layer1, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(e_c_3_layer2, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(e_c_3_layer3, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(e_c_3_layer4, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(e_c_3_layer5, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(e_c_3_layer6, 0.0012f, 3.0f, -68.5f, 10.0f); // FS

	//
	int e_c_5_layer1=sim->createGroup("ec5_1", ceil(neuronsPerGroup*.233), EXCITATORY_NEURON);
	int e_c_5_layer2=sim->createGroup("ec5_2", ceil(neuronsPerGroup*.133), EXCITATORY_NEURON);
	int e_c_5_layer3=sim->createGroup("ec5_3", ceil(neuronsPerGroup*.2), EXCITATORY_NEURON);
	int e_c_5_layer4=sim->createGroup("ec5_4", ceil(neuronsPerGroup*.1), EXCITATORY_NEURON);
	int e_c_5_layer5=sim->createGroup("ec5_5", ceil(neuronsPerGroup*.067), EXCITATORY_NEURON);
	int e_c_5_layer6=sim->createGroup("ec5_6", ceil(neuronsPerGroup*.267), EXCITATORY_NEURON);
	sim->setNeuronParameters(e_c_5_layer1, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(e_c_5_layer2, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(e_c_5_layer3, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(e_c_5_layer4, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(e_c_5_layer5, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(e_c_5_layer6, 0.0012f, 3.0f, -68.5f, 10.0f); // FS

	//
	int c_a_1_layer1=sim->createGroup("ca1_1", ceil(neuronsPerGroup*.233), EXCITATORY_NEURON);
	int c_a_1_layer2=sim->createGroup("ca1_2", ceil(neuronsPerGroup*.133), EXCITATORY_NEURON);
	int c_a_1_layer3=sim->createGroup("ca1_3", ceil(neuronsPerGroup*.2), EXCITATORY_NEURON);
	int c_a_1_layer4=sim->createGroup("ca1_4", ceil(neuronsPerGroup*.1), EXCITATORY_NEURON);
	int c_a_1_layer5=sim->createGroup("ca1_5", ceil(neuronsPerGroup*.067), EXCITATORY_NEURON);
	int c_a_1_layer6=sim->createGroup("ca1_6", ceil(neuronsPerGroup*.267), EXCITATORY_NEURON);
	sim->setNeuronParameters(c_a_1_layer1, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(c_a_1_layer2, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(c_a_1_layer3, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(c_a_1_layer4, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(c_a_1_layer5, 0.0012f, 3.0f, -68.5f, 10.0f); // FS
	sim->setNeuronParameters(c_a_1_layer6, 0.0012f, 3.0f, -68.5f, 10.0f); // FS

	// random connection with 10% probability
	//int c0=sim->connect(spike_gen, e_c_3_layer, "random", RangeWeight(0.005f), 0.1f, RangeDelay(1,10));
	int c0=sim->connect(spike_gen, e_c_3_layer, "random", RangeWeight(0.005f), 1.0);
	//int c1=sim->connect(e_c_3_layer, e_c_5_layer, "random", RangeWeight(0.001f), 0.1f, RangeDelay(1,10));
	int c1=sim->connect(e_c_3_layer, e_c_5_layer, "random", RangeWeight(0.005f), 0.0);
	//int c2=sim->connect(e_c_5_layer, c_a_1_layer, "random", RangeWeight(0.003f), 0.1f, RangeDelay(1,10));
	int c2=sim->connect(e_c_5_layer, c_a_1_layer, "random", RangeWeight(0.005f), 0.0);

	sim->setConductances(false);


	// ---------------- SETUP STATE -------------------

	sim->setupNetwork();

	sim->setExternalCurrent(e_c_3_layer, -190.0);
	sim->setExternalCurrent(e_c_5_layer, -190.0);
	sim->setExternalCurrent(c_a_1_layer, -190.0);

	SpikeMonitor* SpikeMonInput  = sim->setSpikeMonitor(spike_gen,"DEFAULT");

	// accept firing rates within this range of target firing
	double target_firing_e_c_3 = 3.6;//27.4;	// target firing rate for gec3
	double target_firing_e_c_5 = 2.6;//42.8;	// target firing rate for gec5
	double target_firing_c_a_1 = 1.2;//42.8;	// target firing rate for gca1

	// algorithm will terminate when at least one of the termination conditions is reached
	double errorMarginHz = 0.015;	// error margin
	int maxIter = 5;//100;				// max number of iterations

	// set up weight tuning from input -> EC3
	SimpleWeightTuner swt_sg_to_ec3(sim, errorMarginHz, maxIter);
	swt_sg_to_ec3.setConnectionToTune(c0, 0.0); // start at 0
	swt_sg_to_ec3.setTargetFiringRate(e_c_3_layer, target_firing_e_c_3);

	// set up weight tuning from EC3 -> EC5
	SimpleWeightTuner swt_ec3_to_ec5(sim, errorMarginHz, maxIter);
	swt_ec3_to_ec5.setConnectionToTune(c1, 0.0); // start at 0
	swt_ec3_to_ec5.setTargetFiringRate(e_c_5_layer, target_firing_e_c_5);

	// set up weight tuning from EC5 -> CA1
	SimpleWeightTuner swt_ec5_to_ca1(sim, errorMarginHz, maxIter);
	swt_ec5_to_ca1.setConnectionToTune(c2, 0.0); // start at 0
	swt_ec5_to_ca1.setTargetFiringRate(c_a_1_layer, target_firing_c_a_1);


	// ---------------- RUN STATE -------------------

	/*printf("\nMemory module synaptic strength tuning\n");
	printf("- Tune weights from spike generator to EC3\n");
	while (!swt_sg_to_ec3.done()) {
		swt_sg_to_ec3.iterate();
	}

	printf("- Tune weights from EC3 to EC5\n");
		while (!swt_ec3_to_ec5.done()) {
			swt_ec3_to_ec5.iterate();
	}

	printf("- Tune weights from EC5 to CA1\n");
		while (!swt_ec5_to_ca1.done()) {
			swt_ec5_to_ca1.iterate();
	}*/

	printf("\n- Verify result (gec3=%.4fHz, gec5=%.4fHz, gac1=%.4fHz, +/- %.4fHz)\n",
			target_firing_e_c_3, target_firing_e_c_5, target_firing_c_a_1, errorMarginHz);
	sim->runNetwork(10,0);

	delete sim;
	return 0;
}
