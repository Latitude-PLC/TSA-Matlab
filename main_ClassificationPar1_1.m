function main_ClassificationPar1_1(task,ntasks)
% Matlab code for Continous Classification (same as main_ClassPar with task
% and ntasks as input variables)
%
% Revisions: $ Date: 11/25/2015 $ Copyright: Zhe Zhu
% Version 1.1: Combine classification and change detection output (07/05/2015)
% Version 1.0: Fast classification for each line (01/10/2015)
%
% Get variable and path inputs
v_input = ccdc_Inputs;
% 1. number of nrows processed
nrows = v_input.ijdim(1);
% 2. number of pixels procesed per line
ncols = v_input.ijdim(2);
% 3. total number of CPU used
% ntasks = 2; 
% 4. CPU identification number for parallel computing
% task = 1;
% 5. locations of data
dir_l = v_input.l_dir;

%% Constants: 
% number of trees
ntrees = 500;
% version of CCDC
ccdc_v = 1.1;

% log for CCDC Change paramters and versions
% report only for the first task
if task == 1
    fileID = fopen('CCDC_Classification_log.txt','w');
    % write location of image stack
    fprintf(fileID,'Image location = %s\r\n',dir_l);
    % write number of images used
    fprintf(fileID,'Number of trees = %d\r\n',ntrees);
    % CCDC Version
    fprintf(fileID,'CCDC Classification Version = %.2f\r\n',ccdc_v);
    % updates
    fprintf(fileID,'******************************************************************************************************\r\n');
    fprintf(fileID,'Revisions: $ Date: 11/20/2015 $ Copyright: Zhe Zhu\r\n');
    fprintf(fileID,'Version 1.1   Combine classification and change detection output (07/05/2015)\r\n');
    fprintf(fileID,'Version 1.0   Fast classification for each line (01/10/2015)\r\n');
    fprintf(fileID,'******************************************************************************************************\r\n');
    fclose(fileID);
end

% load in RF models 
load modelRF % change name if neccessary

% make TSFitMap folder for storing coefficients
n_result = v_input.name_rst;% 'TSFitMap';
if isempty(dir(n_result))
    mkdir(n_result);
end

% prepare the irows for task for ALL rows
irows = zeros(1,1);
i = 0;
while task + ntasks*i <= nrows
   irows(i+1) = task + ntasks*i;
   i = i+1;
end

for i = 1:length(irows)   
    try
        % if successfully loaded => skip
        load([dir_l,'/',n_result,'/','record_change',sprintf('%d',irows(i)),'.mat']);
    catch me
        % not exist or corrupt
        fprintf('Missing the %dth row!\n',irows(i));
        continue
    end
    
    fprintf('Processing the %dth row\n',irows(i));
    % Continous Classfication Done for a Line of timeseries pixels
    Class_Line1_1(dir_l,n_result,irows(i),modelRF,v_input.num_c,v_input.nbands,ntrees);
end
fprintf('Done!\n');
exit;
end