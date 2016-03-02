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
//#include <array>
//#include <string>
//#include <algorithm>

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

struct create_syn_weight_variables
{
	std::string syn_type;
	int groups_in_layer;
	double inital_firing[];
	double group_sizes[1000];
	double neuronsPerGroup;
	double target_firing_rates[]; // accept firing rates within this range of target firing
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

			/*std::cout<<"\n new connection:\n";
			std::cout<<syn_connections_formed;
			//std::cout<<"\n";
			std::cout<<"\ni:\t";std::cout<<i;std::cout<<"\tinput:\t";std::cout<<SSTR(input_layer[i]);std::cout<<"\toutput:\t";std::cout<<SSTR(output_layer[i]);
			std::cout<<"\n remaining_connections:\n";
			std::cout<<remaining_connections;
			std::cout<<"\n connection_probability: \n";
			std::cout<<connection_probability;
			std::cout<<"\n full_syn_weight: \n";
			std::cout<<adjusted_syn_weight;*/

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

double create_syn_weights(create_syn_weight_variables syn_weight_variables, int i) {
	/*
	 * Synapse weights are generated based on a model fitted to example data
	 * using this tool: http://www.xuru.org/rt/MLR.asp#CopyPaste
	 */
	double synapse_weight;
	double x_1;
	double x_2;
	double x_3;
	double x_4;

	if (syn_weight_variables.syn_type == "ec3_to_ec5") {
		x_1 = 3.912113944*.000001;
		x_2 = -4.16315869*.0001;
		x_3 = 2.623349094*.01;
		x_4 = 2.139377967*.001;
	}
	else if (syn_weight_variables.syn_type == "ec3_to_ec5") {
		x_1 = 3.033292035*.00001;
		x_2 = -2.514693722*.001;
		x_3 = 7.446468689*.01;
		x_4 = -5.280723543*.01;
	}
	//for (int i = 0; i < syn_weight_variables.groups_in_layer; i++) {
		synapse_weight = x_1*syn_weight_variables.inital_firing[i]+x_2*
		(syn_weight_variables.neuronsPerGroup*syn_weight_variables.group_sizes[i])+
		x_3*syn_weight_variables.target_firing_rates[i]+x_4;
	//}
	return synapse_weight;
}

int main(int argc, const char* argv[]) {
	/*
	 * the _conn arrays set the synapse connection quantities
	 */
	double *synapse_weights;
	double sg_to_ec3_initial_firing[] = {5968, 3352, 5002, 2600, 3203, 11137};
	double ec3_to_ec5_initial_firing[] = {8545, 6835, 10002, 1252, 1122, 4154};
	double ec5_to_ca1_initial_firing[] = {56752, 33250, 15803, 8151, 5610, 32574};
	double sg_to_ec3_target_firing[] = {1.4917, 2.2081, 2.2081, 0.6152, 0.3024, 0.3024};
	double ec3_to_ec5_target_firing[] = {1.4917, 2.2081, 2.2081, 0.6152, 0.3024, 0.3024};
	double ec5_to_ca1_target_firing[] = {6.8895, 4.6546, 1.6016, 5.7480, 5.7480, 7.6722};

	// ---------------- CONFIG STATE -------------------
	sim = new CARLsim("MemModGPU", GPU_MODE, USER, 0, 42);

	create_layers_variables sg_layer;
	create_syn_variables sg_to_ec3_synapes;
	create_syn_variables ec3_to_ec5_synapes;
	create_syn_variables ec5_to_ca1_synapes;
	create_syn_weight_variables syn_weight_variables;

	// SpikeGenerator to help feed input to ec3 to setup the simulated layer.
	PeriodicSpikeGenerator PSG_for_ec3(10.0f);//(00.1f);//(50.0f);
	int psg_input = sim->createSpikeGeneratorGroup("psg1",
			sg_layer.neuronsPerGroup, EXCITATORY_NEURON);
	sim->setSpikeGenerator(psg_input, &PSG_for_ec3);

	create_layers_variables e_c_3_layer;
	e_c_3_layer = create_layers(e_c_3_layer);
	create_layers_variables e_c_5_layer;
	e_c_5_layer = create_layers(e_c_5_layer);
	create_layers_variables c_a_1_layer;
	c_a_1_layer = create_layers(c_a_1_layer);

	double sg_to_ec3_conn[sg_to_ec3_synapes.groups_in_layer] = {0.000075, 0.0001, 0.000061, 0.000078, 0.001, 0.0008};
	for (int i = 0; i < sg_to_ec3_synapes.connections_per_group; i++) {sg_to_ec3_synapes.connections_to_form[i]=sg_to_ec3_conn[i];};
	sg_to_ec3_synapes = create_syn(sg_layer.layers, e_c_3_layer.layers, sg_to_ec3_synapes);

	syn_weight_variables.groups_in_layer = sg_to_ec3_synapes.groups_in_layer;
	syn_weight_variables.neuronsPerGroup = sg_layer.neuronsPerGroup;
	syn_weight_variables.syn_type ="sg_to_ec3";
	for (int i = 0; i < sg_layer.groups_in_layer; i++) {syn_weight_variables.group_sizes[i]=sg_layer.group_sizes[i];};
	for (int i = 0; i < sg_layer.groups_in_layer; i++) {syn_weight_variables.inital_firing[i]=ec3_to_ec5_initial_firing[i];};
	for (int i = 0; i < sg_layer.groups_in_layer; i++) {syn_weight_variables.target_firing_rates[i]=ec3_to_ec5_target_firing[i];};
	//synapse_weights = create_syn_weights(syn_weight_variables);

	for (int i = 0; i < ec3_to_ec5_synapes.connections_per_group; i++) {ec3_to_ec5_synapes.connections_to_form[i]=
			create_syn_weights(syn_weight_variables, i);}
	ec3_to_ec5_synapes = create_syn(e_c_3_layer.layers, e_c_5_layer.layers, ec3_to_ec5_synapes);

	for (int i = 0; i < sg_layer.groups_in_layer; i++) {syn_weight_variables.inital_firing[i]=ec5_to_ca1_initial_firing[i];};
	for (int i = 0; i < sg_layer.groups_in_layer; i++) {syn_weight_variables.target_firing_rates[i]=ec5_to_ca1_target_firing[i];};
	//synapse_weights = create_syn_weights(syn_weight_variables);

	for (int i = 0; i < ec5_to_ca1_synapes.connections_per_group; i++) {ec5_to_ca1_synapes.connections_to_form[i]=
			create_syn_weights(syn_weight_variables, i);}
	ec5_to_ca1_synapes = create_syn(e_c_5_layer.layers, c_a_1_layer.layers, ec5_to_ca1_synapes);

	sim->setConductances(false);

	// ---------------- SETUP STATE -------------------

	sim->setupNetwork();

	create_external_current(e_c_3_layer.layers, -172.0, 6);//-172.0, 6);//-171.196, 6);
	create_external_current(e_c_5_layer.layers, -180.0, 6);
	create_external_current(c_a_1_layer.layers, -180.0, 6);

	create_spike_monitors(e_c_3_layer.layers, 6);
	create_spike_monitors(e_c_5_layer.layers, 6);
	create_spike_monitors(c_a_1_layer.layers, 6);

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
	sim->runNetwork(20,0);

	/*sim->biasWeights(c1, step_size, true);*/

	//sim->setWeight(c1, e_c_3_layer1, e_c_3_layer5, 10000.0, true);
	//sim->scaleWeights(c1, 1000.0, true);
	//sim->runNetwork(2,0);
	//sim->runNetwork(1,500);
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

	///// output expected firing ////////
	/*std::cout<<"spike mon output:\t";
	//std::cout<<
	//std::cout<<sim->getSpikeMonitor(0)->getAllFiringRates();//getNeuronNumSpikes(0);
	double thisRate = sim->getSpikeMonitor(e_c_5_layer.layers[0])->getNeuronNumSpikes(0);//getPopMeanFiringRate();
	std::cout<<thisRate;
	std::cout<<"\noutput end";

	double e_c_3_exper_firing_rates[6] = {284.17424, 284.17424, 284.17424, 284.17424,  578.0694,  578.0694};
	//double e_c_5_layer_exper_sizes[6] = {7.0, 4.0, 6.0, 3.0, 2.0, 8.0};
	double e_c_5_exper_firing_rates[6] = {423.9286, 627.5, 627.5, 174.8462, 174.8462, 174.8462};
	//double c_a_1_layer_exper_sizes[6] = {7.0, 4.0, 6.0, 3.0, 2.0, 8.0};
	double c_a_1_exper_firing_rates[6] = {2920.7727, 2920.7727, 1005.0201, 1005.0201, 1005.0201, 1341.4554};
	double orig_exper_neuron_size = 30; //neurons in orig groups
	double time_conversion = 0.01; // 200k/2000k ms simulated from orig exper
	for (int i = 0; i < 6; i++) {
		std::cout<<"\ne_c_3_layer group:";
		std::cout<<i;
		std::cout<<"\tsize:\t";
		std::cout<<e_c_3_layer.group_sizes[i]*e_c_3_layer.neuronsPerGroup;
		std::cout<<"\tfiring:\t";
		std::cout<<(e_c_3_layer.neuronsPerGroup/orig_exper_neuron_size)*e_c_3_exper_firing_rates[i]*time_conversion;
		std::cout<<"\ttotal firing:\t";
		std::cout<<(e_c_3_layer.neuronsPerGroup/orig_exper_neuron_size)*e_c_3_exper_firing_rates[i]*time_conversion*e_c_3_layer.group_sizes[i]*e_c_3_layer.neuronsPerGroup;
	}
	for (int i = 0; i < 6; i++) {
		std::cout<<"\ne_c_5_layer group:";
		std::cout<<i;
		std::cout<<"\tsize:\t";
		std::cout<<e_c_5_layer.group_sizes[i]*e_c_5_layer.neuronsPerGroup;
		std::cout<<"\tfiring:\t";
		std::cout<<(e_c_5_layer.neuronsPerGroup/orig_exper_neuron_size)*e_c_5_exper_firing_rates[i]*time_conversion;
		std::cout<<"\ttotal firing:\t";
		std::cout<<(e_c_5_layer.neuronsPerGroup/orig_exper_neuron_size)*e_c_5_exper_firing_rates[i]*time_conversion*e_c_5_layer.group_sizes[i]*e_c_5_layer.neuronsPerGroup;
	}
	for (int i = 0; i < 6; i++) {
		std::cout<<"\nc_a_1_layer group:";
		std::cout<<i;
		std::cout<<"\tsize:\t";
		std::cout<<c_a_1_layer.group_sizes[i]*c_a_1_layer.neuronsPerGroup;
		std::cout<<"\tfiring:\t";
		std::cout<<(c_a_1_layer.neuronsPerGroup/orig_exper_neuron_size)*c_a_1_exper_firing_rates[i]*time_conversion;
		std::cout<<"\ttotal firing:\t";
		std::cout<<(c_a_1_layer.neuronsPerGroup/orig_exper_neuron_size)*c_a_1_exper_firing_rates[i]*time_conversion*c_a_1_layer.group_sizes[i]*c_a_1_layer.neuronsPerGroup;
	}*/

	delete sim;
	return 0;
}
