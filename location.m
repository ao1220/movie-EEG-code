%% load location data
clear;clc;
InputPath = '';
filelist = dir(InputPath);
outputpath = '';
for i = 1:length(filelist)-2
    filepath = strcat(InputPath,'\', filelist(i+2).name);
    temp = dir(filepath);
    filename = [temp(3).name,'.set'];
    EEG = readEEG([filepath,'\',temp(3).name],filename);
    EEG = eeg_checkset( EEG );
    EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp');
    EEG = eeg_checkset( EEG );
    mkdir(outputpath,filelist(i+2).name);
    output = [outputpath,'/',filelist(i+2).name];
    EEG = pop_saveset( EEG, 'filename',filelist(i+2).name,'filepath',output);
end
