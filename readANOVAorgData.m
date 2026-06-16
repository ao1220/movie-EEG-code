function [network1_1, network1_2, network2_1, network2_2, y1_1, y1_2, y2_1, y2_2, chan] = readANOVAorgData(N1, N2, numofchan, numofbands)
% Load network data for ANOVA.

inputpath1_1 = '\movie1_MCI';
inputpath1_2 = '\movie2_MCI';
inputpath2_1 = '\movie1_HC';
inputpath2_2 = '\movie2_HC';

listing1_1 = dir(inputpath1_1);
listing1_2 = dir(inputpath1_2);
listing2_1 = dir(inputpath2_1);
listing2_2 = dir(inputpath2_2);

for i = 1:N1
    disp(i);
    data_dir1_1 = [inputpath1_1, '\', listing1_1(i+2).name];
    data_dir1_2 = [inputpath1_2, '\', listing1_2(i+2).name];
    load(data_dir1_1);
    network1_1(i,:,:,:) = EEG_results.M_zscore_mean;
    load(data_dir1_2);
    network1_2(i,:,:,:) = EEG_results.M_zscore_mean;
end

for i = 1:N2
    disp(i);
    data_dir2_1 = strcat(inputpath2_1, '\', listing2_1(i+2).name);
    data_dir2_2 = strcat(inputpath2_2, '\', listing2_2(i+2).name);
    load(data_dir2_1);
    network2_1(i,:,:,:) = EEG_results.M_zscore_mean;
    load(data_dir2_2);
    network2_2(i,:,:,:) = EEG_results.M_zscore_mean;
end

network1_1(isinf(network1_1)) = 0;
network1_2(isinf(network1_2)) = 0;
network2_1(isinf(network2_1)) = 0;
network2_2(isinf(network2_2)) = 0;

y1_1 = reshape(network1_1, N1, numofchan^2, numofbands);
y1_2 = reshape(network1_2, N1, numofchan^2, numofbands);
y2_1 = reshape(network2_1, N2, numofchan^2, numofbands);
y2_2 = reshape(network2_2, N2, numofchan^2, numofbands);
chan = EEG_results.parameter.chanlocs;
end
