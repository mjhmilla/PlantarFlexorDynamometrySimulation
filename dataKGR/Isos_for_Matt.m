close all 
clear all
clc

path_c3d = 'I:\l01\AG_BIOMECHANIK\Denis Holzer\$$_DATA\10_JP_01\C3D\MVC\';


frequency_analogs = 1000;

%% Open every file in the folder (path_c3d) one after another
files = dir([path_c3d '\*.c3d']);

for j = 1 : size(files)
    trialname = files(j).name;                                               %filename mit .c3d ende
    trial = trialname(1:end-4);                                              %filename ohne filetype
    % OPEN FILE
    fid = fopen([path_c3d '\' trialname]);                                   %macht eine id für die file
    % CLOSE FILE ID
    fclose(fid);                                                            %schliesst file id
        
        
    %% Get c3d Data
    data_c3d.(trial) = btkReadAcquisition([path_c3d '\' trialname]);
    %get analogs
    analogs.(trial) = btkGetAnalogs(data_c3d.(trial));
    timestamp_analogs.(trial) = ((1:length(analogs.(trial).Electric_Potential_GL))-1)/frequency_analogs;
    torque.(trial) = analogs.(trial).Torque_1*(-1000);  
    
    
    TABEL = table(timestamp_analogs.(trial)',torque.(trial));
    
    writetable(TABEL,'H:\Matlab_Tools\$$_Code\KGR\Isos_for_Matt.xlsx','Sheet',trial)
end

save('H:\Matlab_Tools\$$_Code\KGR\Isos_for_Matt','timestamp_analogs','torque')