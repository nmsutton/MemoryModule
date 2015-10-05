function neuron_numbers = neuron_counts(data_file, max_electrode, max_time)
%% Example run: neuron_counts('ec013.156',8,2000000)
%% ans =   15   14   10   13    5    5    4    1
%%

neuron_numbers = [];
for i = [1:max_electrode]
  [T,G,Map,Par]=LoadCluRes(data_file,i);
  clusters = unique(G)';
  neuron_numbers = [neuron_numbers, length(clusters)];
end

%neuron_numbers = sum(neuron_numbers)