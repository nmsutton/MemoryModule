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
#include <sstream>
#include <iostream>
#include <stdlib.h>
//#include <algorithm>

#if defined(WIN32) || defined(WIN64)
	#define _CRT_SECURE_NO_WARNINGS
	#include <periodic_spikegen.h>
	#include <simple_weight_tuner.h>
#endif

/* Method to convert int to str
 * from http://stackoverflow.com/questions/5590381/easiest-way-to-convert-int-to-string-in-c
 */
#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

int layers_formed = 0;
int syn_connections_formed = 0;
CARLsim *sim;

struct create_layers_variables
/*
 * Parameters for creating layers.  Neuron sub groups are sized according to group_sizes[].
 * neuron_parameters[] specify izhikevich neuron parameters.
 * TODO:  Currently hardcoded multiple max parameter numbers as as layer_sim_groups for now.  For the future
 * work on having them dynamically sized.
 */
{
	int layers[1000];
	int layer_sim_groups[1000];
	double group_sizes[1000] = {0.233, 0.133, 0.2, 0.1, 0.067, 0.267};
	double neuron_parameters[4] = {0.0012f, 3.0f, -68.5f, 10.0f};
	double groups_in_layer = 6;
	double neuronsPerGroup = 500;
	CARLsim *sim;
};

struct create_syn_variables
/*
 * Parameters for creating synapses
 */
{
	double fire_rate_ratios[];
	double syn_weights[];
	static const int groups_in_layer = 6;
	double connections_per_group = 8;
	double connections_to_form[groups_in_layer] = {0.1, 10.0, 0.1, 0.1, 0.1, 0.1};
	int syn_connections[10000];
	CARLsim *sim;
};

create_layers_variables create_layers(create_layers_variables layers_variables) {
	/*
	 * Generate a new layer of neurons.
	 */
	int initial_layers_formed = layers_formed;
	int last_new_layer_index = layers_formed + layers_variables.groups_in_layer;
	int normalized_group_index = 0;
	for (int i = layers_formed; i < last_new_layer_index; i++) {
		normalized_group_index = i - initial_layers_formed;
		layers_variables.layers[normalized_group_index]= sim->createGroup(SSTR(i),
				ceil(layers_variables.neuronsPerGroup*layers_variables.group_sizes[normalized_group_index]), EXCITATORY_NEURON);
		sim->setNeuronParameters(layers_variables.layers[normalized_group_index], layers_variables.neuron_parameters[0],
				layers_variables.neuron_parameters[1], layers_variables.neuron_parameters[2], layers_variables.neuron_parameters[3]); // FS
		//std::cout<<layers_formed;
		//std::cout<<"\n";
		layers_formed++;
	}

	return layers_variables;
}

create_syn_variables create_syn(int input_layer[1000], int output_layer[1000], create_syn_variables syn_variables) {
	/*
	 * Generate synaptic connections between layers.
	 * If last section of connections is a fraction then reduce the weight of that section
	 * proportional to the fraction of connections (e.g. create 6.5 out of 8.0 connections, the 7th
	 * section of connections has 0.5/1.0 of the full weight.
	 *
	 * normalized_index = connections index normalized to have the first index at 0
	 * normalized_delay = delay index normalized to have the first index at 1.5
	 */
	//int new_connections_last_index = syn_variables.groups_in_layer + syn_connections_formed;
	int normalized_delay = 0;
	//int initial_syn_connections_formed = syn_connections_formed;
	//int normalized_index = 0;
	double remaining_connections = 0;
	double initial_syn_weight = 10.0f;
	double adjusted_syn_weight = 0.0f;
	double connection_probability = 1.0;
	//bool create_connection = true;

	for (int i = 0; i < syn_variables.groups_in_layer; i++) {
		adjusted_syn_weight = initial_syn_weight;
		connection_probability = 1.0;
		for (int i2 = 0; i2 < syn_variables.connections_per_group; i2++) {
			//normalized_index = i - initial_syn_connections_formed;
			//double times_greater_ratio = int(ceil(syn_variables.fire_rate_ratios[i]));
			normalized_delay = 1.5 + i2;//(new_connections_last_index - i2);

			// Check for last connection section
			remaining_connections = syn_variables.connections_to_form[i] - i2;
			if (remaining_connections < 1.0 & remaining_connections > 0.0) {
				adjusted_syn_weight = initial_syn_weight * remaining_connections;
			}
			else if (remaining_connections <= 0.0) {
				connection_probability = 0.0;
			}

			syn_variables.syn_connections[syn_connections_formed]=sim->connect(input_layer[i],
					output_layer[i], "full", RangeWeight(adjusted_syn_weight), connection_probability, normalized_delay);

			std::cout<<"\n new connection:\n";
			std::cout<<syn_connections_formed;
			//std::cout<<"\n";
			std::cout<<"\ni:\t";std::cout<<i;std::cout<<"\tinput:\t";std::cout<<SSTR(input_layer[i]);std::cout<<"\toutput:\t";std::cout<<SSTR(output_layer[i]);
			std::cout<<"\n remaining_connections:\n";
			std::cout<<remaining_connections;
			std::cout<<"\n connection_probability: \n";
			std::cout<<connection_probability;
			std::cout<<"\n full_syn_weight: \n";
			std::cout<<adjusted_syn_weight;

			syn_connections_formed++;
		}
	}

	return syn_variables;
}

void create_external_current(int layer[1000], int current_value, int groups_to_use) {
	for (int i = 0; i < groups_to_use; i++) {
		sim->setExternalCurrent(layer[i], current_value);
	}
}

void create_spike_monitors(int layer[1000], int groups_to_use) {
	for (int i = 0; i < groups_to_use; i++) {
		SpikeMonitor* SpikeMonInput2  = sim->setSpikeMonitor(layer[i],"DEFAULT");
	}
}

int main(int argc, const char* argv[]) {
	// ---------------- CONFIG STATE -------------------
	sim = new CARLsim("MemModGPU", GPU_MODE, USER, 0, 42);
	int neuronsPerGroup = 500;//500;

	create_layers_variables e_c_3_layer;
	e_c_3_layer = create_layers(e_c_3_layer);
	create_layers_variables e_c_5_layer;
	e_c_5_layer = create_layers(e_c_5_layer);
	create_layers_variables c_a_1_layer;
	c_a_1_layer = create_layers(c_a_1_layer);

	create_syn_variables ec3_to_ca5_synapes;
	//ec3_to_ca5_synapes.connections_to_form = 0.1;//6.5;
	ec3_to_ca5_synapes = create_syn(e_c_3_layer.layers, e_c_5_layer.layers, ec3_to_ca5_synapes);

	create_syn_variables ec5_to_ca1_synapes;
	//ec5_to_ca1_synapes.connections_to_form = 0.1;//4.5;
	ec5_to_ca1_synapes = create_syn(e_c_5_layer.layers, c_a_1_layer.layers, ec5_to_ca1_synapes);

	sim->setConductances(false);

	// ---------------- SETUP STATE -------------------

	sim->setupNetwork();

	create_external_current(e_c_3_layer.layers, -160.0, 6);
	create_external_current(e_c_5_layer.layers, -180.0, 6);
	create_external_current(c_a_1_layer.layers, -180.0, 6);

	create_spike_monitors(e_c_3_layer.layers, 6);
	create_spike_monitors(e_c_5_layer.layers, 3);
	create_spike_monitors(c_a_1_layer.layers, 1);

	// accept firing rates within this range of target firing
	double target_firing_e_c_3_1 = 1.49;//27.4;	// target firing rate for gec3
	double target_firing_e_c_5 = 2.6;//42.8;	// target firing rate for gec5
	double target_firing_c_a_1 = 1.2;//42.8;	// target firing rate for gca1

	// algorithm will terminate when at least one of the termination conditions is reached
	double errorMarginHz = 0.015;	// error margin
	int maxIter = 30;//100;				// max number of iterations
	double startingWeight = 0.05;
	double stepSize = 1.0;

	// set up weight tuning from input -> EC3
/*
	SimpleWeightTuner swt_sg_to_ec3(sim, errorMarginHz, maxIter, stepSize);
	swt_sg_to_ec3.setConnectionToTune(c0, startingWeight, true); // start at 0
	swt_sg_to_ec3.setTargetFiringRate(e_c_3_layer1, target_firing_e_c_3_1);
*/
/*
	// set up weight tuning from EC3 -> EC5
	SimpleWeightTuner swt_ec3_to_ec5(sim, errorMarginHz, maxIter);
	swt_ec3_to_ec5.setConnectionToTune(c1, 0.0); // start at 0
	swt_ec3_to_ec5.setTargetFiringRate(e_c_5_layer, target_firing_e_c_5);

	// set up weight tuning from EC5 -> CA1
	SimpleWeightTuner swt_ec5_to_ca1(sim, errorMarginHz, maxIter);
	swt_ec5_to_ca1.setConnectionToTune(c2, 0.0); // start at 0
	swt_ec5_to_ca1.setTargetFiringRate(c_a_1_layer, target_firing_c_a_1);
*/

	// ---------------- RUN STATE -------------------
/*
	printf("\nMemory module synaptic strength tuning\n");
	printf("- Tune weights from spike generator to EC3\n");
	while (!swt_sg_to_ec3.done()) {
		swt_sg_to_ec3.iterate();
	}
	*/
/*
	printf("- Tune weights from EC3 to EC5\n");
		while (!swt_ec3_to_ec5.done()) {
			swt_ec3_to_ec5.iterate();
	}

	printf("- Tune weights from EC5 to CA1\n");
		while (!swt_ec5_to_ca1.done()) {
			swt_ec5_to_ca1.iterate();
	}*/

	//printf("\n- Verify result (gec3=%.4fHz, gec5=%.4fHz, gac1=%.4fHz, +/- %.4fHz)\n",
	//		target_firing_e_c_3, target_firing_e_c_5, target_firing_c_a_1, errorMarginHz);
	double step_size = 100.0;
	sim->runNetwork(2,0);

	/*sim->biasWeights(c1, step_size, true);*/

	//sim->setWeight(c1, e_c_3_layer1, e_c_3_layer5, 10000.0, true);
	//sim->scaleWeights(c1, 1000.0, true);
	sim->runNetwork(2,0);
	sim->runNetwork(1,500);
	/*sim->setWeight(c0, spike_gen, e_c_3_layer1, 0.0, false);
	sim->scaleWeights(c0, 0.9, false);
	sim->runNetwork(2,0);
	sim->setWeight(c0, spike_gen, e_c_3_layer1, 0.0, false);
	sim->scaleWeights(c0, 2.0, false);
	sim->runNetwork(2,0);
	sim->setWeight(c0, spike_gen, e_c_3_layer1, 3.8, false);
	sim->scaleWeights(c0, 0.8, false);
	sim->runNetwork(2,0);
	sim->setWeight(c0, spike_gen, e_c_3_layer1, 0.0, false);
	sim->scaleWeights(c0, 2.0, false);
	sim->runNetwork(10,0);*/

	delete sim;
	return 0;
}
