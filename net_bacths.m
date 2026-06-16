clear;clc;
inputpath = '';
Output = '';
mkdir(Output);
combs_project_id = '1';
epochLenth = 5;
eventlabel = '[]';
grouplabel = '[]';
method = 'plv';%'psi' 'plv' 'corr' 'ciplv'
proportion = 0;
seleChanns = 'all';
bandLimit = [1,4;4,8;8,12.5;12.5,30;30,40;1,40];
srate = '[]';

filelist = dir(inputpath);
for i = 1:length(filelist)-2
    Input = strcat(inputpath,'\', filelist(i+2).name);
    disp(i)
    wb_pipeline_EEG_calcEEGnetwork(Input,Output,combs_project_id,epochLenth,...
        eventlabel,bandLimit,seleChanns,method,proportion,srate,grouplabel);
end