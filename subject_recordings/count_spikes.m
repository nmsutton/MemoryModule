function spikes = count_spikes(data_file, electrode, max_time)
%% Count the spikes from all neurons in a given dataset 
%% and electrode number
%% max_time = max recording time
%%
%% Example run: count_spikes('ec013.156',5,2000000)
%% ans =  6815

total_spikes = [];

[T,G,Map,Par]=LoadCluRes(data_file,electrode);
clusters = unique(G)';
for i=[clusters]
  clust_times = T(G == i);
  times_filtered = clust_times < max_time;
  total_spikes = [total_spikes, sum(times_filtered)];
end

spikes = sum(total_spikes);
%spikes = (total_spikes);