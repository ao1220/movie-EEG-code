clear;clc;
close all;
InputPath = '';
Output = '';
filelist = dir(InputPath);
mkdir(Output);
combs_project_id = '1';
% para
seleChanns = '[1:31,33:64]';
EOGchanns = '[]';
thre_ODQ = '0';              
passband = '[1,40]';   
PowerFrequency = '50';    
keepUnselectChannsFlag = '0';  
badChannelInterploateFlag = '1'; 
residualArtifactRemovalFlag = '1'; 
MARA_thre = '0.9';        
srate = '1000';
for i = 1:length(filelist)-2
    Input = strcat(InputPath,'\', filelist(i+2).name);
    disp(i)
    wb_pipeline_EEG_prepro(Input,Output,combs_project_id,seleChanns,EOGchanns,...
        thre_ODQ,passband,PowerFrequency,keepUnselectChannsFlag,badChannelInterploateFlag,...
        residualArtifactRemovalFlag,MARA_thre,srate);
end