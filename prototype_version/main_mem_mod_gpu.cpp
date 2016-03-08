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
	double neuron_parameters[4] = {0.0012, 3.0, -65.8, 10.0};
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
	double max_connections_per_group = 20;
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
	int normalized_delay = 0;
	double remaining_connections = 0;
	double initial_syn_weight = 10.0;
	double adjusted_syn_weight = 0.0;
	double connection_probability = 1.0;

	for (int i = 0; i < syn_variables.groups_in_layer; i++) {
		adjusted_syn_weight = initial_syn_weight;
		connection_probability = 1.0;
		for (int i2 = 0; i2 < syn_variables.max_connections_per_group; i2++) {
			normalized_delay = 1.5 + i2;

			// Check for last connection section
			remaining_connections = syn_variables.connections_to_form[i] - i2;
			if (remaining_connections < 1.0 & remaining_connections > 0.0) {
				adjusted_syn_weight = initial_syn_weight * remaining_connections;
			}
			else if (remaining_connections <= 0.0) {
				connection_probability = 0.0;
			}

			if (remaining_connections >= 0.0) {
			syn_variables.syn_connections[syn_connections_formed]=sim->connect(input_layer[i],
					output_layer[i], "full", RangeWeight(adjusted_syn_weight), connection_probability, normalized_delay);

			/*std::cout<<"\n new connection:\n";
			std::cout<<syn_connections_formed;
			//std::cout<<"\n";
			std::cout<<"\ni:\t";std::cout<<i;std::cout<<"\tinput:\t";std::cout<<SSTR(input_layer[i]);std::cout<<"\toutput:\t";std::cout<<SSTR(output_layer[i]);
			std::cout<<"\n remaining_connections:\n";
			std::cout<<remaining_connections;
			std::cout<<"\n connection_probability: \n";
			std::cout<<connection_probability;
			std::cout<<"\n full_syn_weight: \n";
			std::cout<<adjusted_syn_weight;
			std::cout<<"\n output_layer[0]: \n";
			std::cout<<output_layer[0];*/
			}

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

	// y = -3.243224433·10-6 x12 - 8.453549021·10-4 x1 x2 + 8.808208543·10-1 x22 + 3.491527807·10-3 x1 - 6.541997654·10-2 x2 - 5.459357633·10-1
	if (syn_type == "ec3_to_ec5") {
		a = -3.243224433*pow(10, -6);
		b = -8.453549021*pow(10, -4);
		c = 8.808208543*pow(10, -1);
		d = 3.491527807*pow(10, -3);
		e = -6.541997654*pow(10, -2);
		f = -5.459357633*pow(10, -1);
		synapse_weight = a * pow(x_1, 2) + b * (x_1*x_2) + c * pow(x_2, 2) + d * x_1 + e * x_2 + f;
	}
	// y = 3.18151867·10-4 x12 - 1.525851879·10-1 x1 x2 - 11.90618389 x22 + 8.041278944·10-1 x1 + 197.9136262 x2 - 732.4225399
	else if (syn_type == "ec5_to_ca1") {
		a = 3.18151867*pow(10,-4);
		b = -1.525851879*pow(10, -1);
		c = -11.90618389;
		d = 8.041278944*pow(10, -1);
		e = 197.9136262;
		f = -732.4225399;
		synapse_weight = a * pow(x_1, 2) + b * (x_1*x_2) + c * pow(x_2, 2) + d * x_1 + e * x_2 + f;
	}
	cout<<"\n++syn\t";cout<<i;cout<<"\tv\t";cout<<synapse_weight;cout<<"\n";
	cout<<a;cout<<"\t";
	cout<<pow(x_1, 2);cout<<"\t";
	cout<<b;cout<<"\t";cout<<(x_1*x_2);cout<<"\t";cout<<c;cout<<"\t";
	cout<<pow(x_2, 2);cout<<"\t";cout<<d;cout<<"\t";cout<<x_1;cout<<"\t";cout<<e;cout<<"\t";cout<<x_2;cout<<"\t";cout<<f;
	/*cout<<"\n++syn\t";cout<<i;cout<<"\tx1\t";cout<<x_1;cout<<"\tif\t";cout<<init_firing;cout<<"\tx2\t";
	cout<<group_size;cout<<"\tx3\t";
	cout<<x_3*targ_firing;cout<<"\tx4\t";cout<<x_4;cout<<"\tv\t";cout<<synapse_weight;cout<<"\n";*/
	//(syn_weight_variables.neuronsPerGroup*
	return synapse_weight;
}

int main(int argc, const char* argv[]) {
	/*
	 * the _conn arrays set the synapse connection quantities
	 */
	std::cout.precision(17);
	double ec3_to_ec5_initial_firing[] = {351.0, 201.0, 300.0, 150.0, 204.0, 804.0};
	double ec5_to_ca1_initial_firing[] = {468.0, 402.0, 600.0, 100.0, 34.0, 268.0};
	double ec3_to_ec5_target_firing[] = {1.4918, 2.2082, 2.2082, 0.6153, 0.3025, 0.3025};
	double ec5_to_ca1_target_firing[] = {6.8898, 4.6546, 1.6016, 5.7480, 5.7480, 7.6722};

	// ---------------- CONFIG STATE -------------------
	sim = new CARLsim("MemModGPU", GPU_MODE, USER, 0, 42);

	create_layers_variables sg_layer;
	create_syn_variables sg_to_ec3_synapes;
	create_syn_variables ec3_to_ec5_synapes;
	create_syn_variables ec5_to_ca1_synapes;

	// SpikeGenerator to help feed input to ec3 to setup the simulated layer.
	PeriodicSpikeGenerator PSG_for_ec3(0.05);
	int psg_input = sim->createSpikeGeneratorGroup("psg1",
			sg_layer.neuronsPerGroup, EXCITATORY_NEURON);
	sim->setSpikeGenerator(psg_input, &PSG_for_ec3);

	create_layers_variables e_c_3_layer;
	e_c_3_layer = create_layers(e_c_3_layer);
	create_layers_variables e_c_5_layer;
	e_c_5_layer = create_layers(e_c_5_layer);
	create_layers_variables c_a_1_layer;
	c_a_1_layer = create_layers(c_a_1_layer);

	// synapses for sg_to_ec3
	//double sg_to_ec3_conn[sg_to_ec3_synapes.groups_in_layer] = {0.000206975, 0.00020585, 0.0002059, 0.0002059, 0.00020710845, 0.0002071121049};
	double sg_to_ec3_conn[sg_to_ec3_synapes.groups_in_layer] = {2.0, 2.4, 2.6, 2.8, 5.0, 5.0};
	for (int i = 0; i < sg_to_ec3_synapes.groups_in_layer; i++) {sg_to_ec3_synapes.connections_to_form[i]=sg_to_ec3_conn[i];};
	sg_to_ec3_synapes = create_syn(sg_layer.layers, e_c_3_layer.layers, sg_to_ec3_synapes);

	// synapses for ec3_to_ec5
	//double ec3_to_ec5_conn[ec3_to_ec5_synapes.groups_in_layer] = {1.7, 3.8, 3.8, 0.12, 0.04, 0.02};
	for (int i = 0; i < ec3_to_ec5_synapes.groups_in_layer; i++) {ec3_to_ec5_synapes.connections_to_form[i]=//ec3_to_ec5_conn[i];};
			create_syn_weights("ec3_to_ec5", sg_to_ec3_synapes.groups_in_layer, ec3_to_ec5_initial_firing[i], sg_layer.group_sizes[i],
								sg_layer.neuronsPerGroup, ec3_to_ec5_target_firing[i], i);};
	ec3_to_ec5_synapes = create_syn(e_c_3_layer.layers, e_c_5_layer.layers, ec3_to_ec5_synapes);

	// synapses for ec5_to_ca1
	//double ec5_to_ca1_conn[ec5_to_ca1_synapes.groups_in_layer] = {20.0, 20.0, 4.4, 7.7, 9.7, 9.8};
	for (int i = 0; i < ec5_to_ca1_synapes.groups_in_layer; i++) {ec5_to_ca1_synapes.connections_to_form[i]=//ec5_to_ca1_conn[i];};
			create_syn_weights("ec5_to_ca1", ec5_to_ca1_synapes.groups_in_layer, ec5_to_ca1_initial_firing[i], sg_layer.group_sizes[i],
								sg_layer.neuronsPerGroup, ec5_to_ca1_target_firing[i], i);};
	ec5_to_ca1_synapes = create_syn(e_c_5_layer.layers, c_a_1_layer.layers, ec5_to_ca1_synapes);

	sim->setConductances(false);

	// ---------------- SETUP STATE -------------------

	sim->setupNetwork();

	create_external_current(e_c_3_layer.layers, -190.0, 6);
	create_external_current(e_c_5_layer.layers, -180.0, 6);
	create_external_current(c_a_1_layer.layers, -180.0, 6);

	create_spike_monitors(e_c_3_layer.layers, 6);
	create_spike_monitors(e_c_5_layer.layers, 6);
	create_spike_monitors(c_a_1_layer.layers, 6);

	sim->runNetwork(20,0);

	delete sim;
	return 0;
}
