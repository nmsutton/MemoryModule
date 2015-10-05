function spikes_combine = elect_spikes(data_file, max_electrode, max_time)
%% Return spike counts for all electrodes in a dataset

spikes_combine = [];
for i = [1:max_electrode]
  spikes = count_spikes(data_file,i,max_time);
  spikes_combine = [spikes_combine, spikes];
end