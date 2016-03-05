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

using namespace std;

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
int spike_monitors_formed = 0;
SpikeMonitor* spike_monitors[1000];
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
	double connections_to_form[groups_in_layer];
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
	double normalized_delay = 0;
	double remaining_connections = 0;
	double initial_syn_weight = 100.0f;//10.0f;
	double adjusted_syn_weight = 0.0f;
	double connection_probability = 1.0;

	for (int i = 0; i < syn_variables.groups_in_layer; i++) {
		adjusted_syn_weight = initial_syn_weight;
		connection_probability = 1.0;
		for (double i2 = 0; i2 < syn_variables.connections_per_group; i2++) {
			normalized_delay = 1.5 + i2;

			// Check for last connection section
			remaining_connections = syn_variables.connections_to_form[i] - i2;
			if  (remaining_connections >= 1.0) {
				//adjusted_syn_weight = 1.0;
			}
			else if(remaining_connections < 1.0 & remaining_connections > 0.0) {
				adjusted_syn_weight = initial_syn_weight * remaining_connections;
			}
			else if (remaining_connections <= 0.0) {
				connection_probability = 0.0;
			}

			syn_variables.syn_connections[syn_connections_formed]=sim->connect(input_layer[i],
					output_layer[i], "full", RangeWeight(adjusted_syn_weight), connection_probability, normalized_delay);

			syn_connections_formed++;
		}
	}

	return syn_variables;
}

void create_external_current(int layer[1000], double current_value, int groups_to_use) {
	for (int i = 0; i < groups_to_use; i++) {
		sim->setExternalCurrent(layer[i], current_value);
	}
}

void create_spike_monitors(int layer[1000], int groups_to_use) {
	for (int i = 0; i < groups_to_use; i++) {
		spike_monitors[spike_monitors_formed] = sim->setSpikeMonitor(layer[i],"DEFAULT");
		spike_monitors_formed++;
	}
}

double create_syn_weights(string syn_type, int group_number, double init_firing, double group_size,
		double tot_neurons, double targ_firing, int i) {
	/*
	 * Synapse weights are generated based on a model fitted to example data
	 * using this tool: http://www.xuru.org/rt/MPR.asp
	 */
	double synapse_weight;
	double a;
	double b;
	double c;
	double d;
	double e;
	double f;
	double x_1 = init_firing;
	double x_2 = targ_firing;

	// y = -1.179096929·10-9 x12 - 8.597476548·10-6 x1 x2 + 5.91313401·10-3 x22 + 1.922300225·10-5 x1 + 3.464477757·10-2 x2 - 4.90256488·10-2
	if (syn_type == "ec3_to_ec5") {
		a = -1.179096929*pow(10, -9);
		b = -8.597476548*pow(10, -6);
		c = 5.91313401*pow(10, -3);
		d = 1.922300225*pow(10, -5);
		e = 3.464477757*pow(10, -2);
		f = -4.90256488*pow(10, -2);
		synapse_weight = a * pow(x_1, 2) + b * (x_1*x_2) + c * pow(x_2, 2) + d * x_1 + e * x_2 + f;
	}
	// y = 1.880657863·10-8 x12 - 3.280482222·10-4 x1 x2 - 3.985943993·10-1 x22 + 1.879435901·10-3 x1 + 6.690053854 x2 - 25.09180372
	else if (syn_type == "ec5_to_ca1") {
		a = 1.880657863*pow(10,-8);
		b = -3.280482222*pow(10, -4);
		c = -3.985943993*pow(10, -1);
		d = 1.879435901*pow(10, -3);
		e = 6.690053854;
		f = -25.09180372;
		synapse_weight = a * pow(x_1, 2) + b * (x_1*x_2) + c * pow(x_2, 2) + d * x_1 + e * x_2 + f;
	}

	return synapse_weight;
}

/*int * create_spike_generators(int neuronsPerGroup, int groups_in_layer, float sg_vals[], string sg_names[]) {
	//
	// SpikeGenerator to generate ec3's input to setup the first input to the simulation.
	// TODO: Make the creations here more automated instead of the repeated lines.
	//
	int sg_ids[groups_in_layer];
	PeriodicSpikeGenerator PSG_for_ec3_1(20.0f);
	PeriodicSpikeGenerator PSG_for_ec3_2(20.0f);
	PeriodicSpikeGenerator PSG_for_ec3_3(20.0f);
	PeriodicSpikeGenerator PSG_for_ec3_4(20.0f);
	PeriodicSpikeGenerator PSG_for_ec3_5(20.0f);
	PeriodicSpikeGenerator PSG_for_ec3_6(20.0f);
	//for (int i = 0; i < groups_in_layer; i++) {
	//	sg_ids[i] = sim->createSpikeGeneratorGroup(sg_names[i], neuronsPerGroup, EXCITATORY_NEURON);
	//}
	sg_ids[0] = sim->createSpikeGeneratorGroup("test0", neuronsPerGroup, EXCITATORY_NEURON);
	sg_ids[1] = sim->createSpikeGeneratorGroup("test1", neuronsPerGroup, EXCITATORY_NEURON);
	sg_ids[2] = sim->createSpikeGeneratorGroup("test2", neuronsPerGroup, EXCITATORY_NEURON);
	sg_ids[3] = sim->createSpikeGeneratorGroup("test3", neuronsPerGroup, EXCITATORY_NEURON);
	sg_ids[4] = sim->createSpikeGeneratorGroup("test4", neuronsPerGroup, EXCITATORY_NEURON);
	sg_ids[5] = sim->createSpikeGeneratorGroup("test5", neuronsPerGroup, EXCITATORY_NEURON);
	sim->setSpikeGenerator(sg_ids[0], &PSG_for_ec3_1);
	sim->setSpikeGenerator(sg_ids[1], &PSG_for_ec3_2);
	sim->setSpikeGenerator(sg_ids[2], &PSG_for_ec3_3);
	sim->setSpikeGenerator(sg_ids[3], &PSG_for_ec3_4);
	sim->setSpikeGenerator(sg_ids[4], &PSG_for_ec3_5);
	sim->setSpikeGenerator(sg_ids[5], &PSG_for_ec3_6);

	return sg_ids;
}*/

int main(int argc, const char* argv[]) {
	/*
	 * the _conn arrays set the synapse connection quantities
	 */
	std::cout.precision(17);
	double ec3_to_ec5_initial_firing[] = {5968.0, 3352.0, 5002.0, 2600.0, 3203.0, 11137.0};
	double ec5_to_ca1_initial_firing[] = {8545.0, 6835.0, 10002.0, 1252.0, 1122.0, 4154.0};
	double ec3_to_ec5_target_firing[] = {1.4918, 2.2082, 2.2082, 0.6153, 0.3025, 0.3025};
	double ec5_to_ca1_target_firing[] = {6.8898, 4.6546, 1.6016, 5.7480, 5.7480, 7.6722};
	//create_layers_variables sg_layer;
	//create_syn_variables sg_to_ec3_synapes;
	create_layers_variables e_c_5_layer;
	create_layers_variables c_a_1_layer;
	create_syn_variables ec3_to_ec5_synapes;
	create_syn_variables ec5_to_ca1_synapes;
	double ec3[ec3_to_ec5_synapes.groups_in_layer];
	double sg_vals[] = {0.0315, 1894.4949/20000, 1420.8712/20000, 710.4356/20000, 963.4490/20000, 3853.7960/20000};
	string sg_names[] = {"ec3_sg_1", "ec3_sg_2", "ec3_sg_3", "ec3_sg_4", "ec3_sg_5", "ec3_sg_6"};
	int *e_c_3_layer;

	// ---------------- CONFIG STATE -------------------
	sim = new CARLsim("MemModGPU", GPU_MODE, USER, 0, 42);

	//e_c_3_layer = create_spike_generators(e_c_5_layer.neuronsPerGroup, ec3_to_ec5_synapes.groups_in_layer, sg_vals, sg_names);
	int groups_in_layer = ec3_to_ec5_synapes.groups_in_layer;
	int neuronsPerGroup = e_c_5_layer.neuronsPerGroup;
	int sg_ids[groups_in_layer];
	PeriodicSpikeGenerator PSG_for_ec3_1(sg_vals[0]);
	PeriodicSpikeGenerator PSG_for_ec3_2(sg_vals[1]);
	PeriodicSpikeGenerator PSG_for_ec3_3(sg_vals[2]);
	PeriodicSpikeGenerator PSG_for_ec3_4(sg_vals[3]);
	PeriodicSpikeGenerator PSG_for_ec3_5(sg_vals[4]);
	PeriodicSpikeGenerator PSG_for_ec3_6(sg_vals[5]);
	//for (int i = 0; i < groups_in_layer; i++) {
		sg_ids[0] = sim->createSpikeGeneratorGroup(sg_names[0], neuronsPerGroup, EXCITATORY_NEURON);
	//}
	sim->setSpikeGenerator(sg_ids[0], &PSG_for_ec3_1);
	//sim->setSpikeRate(sg_ids[0], )
	/*sim->setSpikeGenerator(sg_ids[1], &PSG_for_ec3_2);
	sim->setSpikeGenerator(sg_ids[2], &PSG_for_ec3_3);
	sim->setSpikeGenerator(sg_ids[3], &PSG_for_ec3_4);
	sim->setSpikeGenerator(sg_ids[4], &PSG_for_ec3_5);
	sim->setSpikeGenerator(sg_ids[5], &PSG_for_ec3_6);*/

	//create_layers_variables e_c_3_layer;
	//e_c_3_layer = create_layers(e_c_3_layer);
	//create_layers_variables e_c_5_layer;
	e_c_5_layer = create_layers(e_c_5_layer);
	//create_layers_variables c_a_1_layer;
	c_a_1_layer = create_layers(c_a_1_layer);

	// synapses for sg_to_ec3
	/*double sg_to_ec3_conn[sg_to_ec3_synapes.groups_in_layer] = {0.00002874007, 0.00002, 0.000028, 0.000008, 0.00002, 0.00001};
	for (int i = 0; i < sg_to_ec3_synapes.groups_in_layer; i++) {sg_to_ec3_synapes.connections_to_form[i]=sg_to_ec3_conn[i];};
	sg_to_ec3_synapes = create_syn(sg_layer.layers, e_c_3_layer.layers, sg_to_ec3_synapes);*/

	// synapses for ec3_to_ec5
	//double ec3_to_ec5_conn[ec3_to_ec5_synapes.groups_in_layer] = {0.012, 0.04386, 0.028, 0.002785, 0.00314, 0.000871};
	//double ec3_to_ec5_conn[ec3_to_ec5_synapes.groups_in_layer] = {0.018, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001};
	double ec3_to_ec5_conn[ec3_to_ec5_synapes.groups_in_layer] = {8.0, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001};
	for (int i = 0; i < ec3_to_ec5_synapes.groups_in_layer; i++) {ec3_to_ec5_synapes.connections_to_form[i]=ec3_to_ec5_conn[i];};
			/*create_syn_weights("ec3_to_ec5", 1, ec3_to_ec5_initial_firing[i], e_c_5_layer.group_sizes[i],
					e_c_5_layer.neuronsPerGroup, ec3_to_ec5_target_firing[i], i);};*/
	ec3_to_ec5_synapes = create_syn(sg_ids, e_c_5_layer.layers, ec3_to_ec5_synapes);

	// synapses for ec5_to_ca1
	double ec5_to_ca1_conn[ec5_to_ca1_synapes.groups_in_layer] = {1.0, 0.7, 0.025, 0.215, 0.21, 0.45};
	for (int i = 0; i < ec5_to_ca1_synapes.groups_in_layer; i++) {ec5_to_ca1_synapes.connections_to_form[i]=ec5_to_ca1_conn[i];};
			/*create_syn_weights("ec5_to_ca1", ec5_to_ca1_synapes.groups_in_layer, ec5_to_ca1_initial_firing[i], sg_layer.group_sizes[i],
								sg_layer.neuronsPerGroup, ec5_to_ca1_target_firing[i], i);};*/
	ec5_to_ca1_synapes = create_syn(e_c_5_layer.layers, c_a_1_layer.layers, ec5_to_ca1_synapes);

	sim->setConductances(false);

	// ---------------- SETUP STATE -------------------

	sim->setupNetwork();

	//create_external_current(e_c_3_layer.layers, -172.0, 6);
	create_external_current(e_c_5_layer.layers, -180.0, 6);
	create_external_current(c_a_1_layer.layers, -180.0, 6);

	//create_spike_monitors(sg_layer.layers, 6);
	SpikeMonitor* SpikeMonInput  = sim->setSpikeMonitor(sg_ids[0],"DEFAULT");
	//create_spike_monitors(e_c_3_layer.layers, 6);
	create_spike_monitors(e_c_5_layer.layers, 1);
	create_spike_monitors(c_a_1_layer.layers, 1);

	sim->runNetwork(200,0);

	delete sim;
	return 0;
}
