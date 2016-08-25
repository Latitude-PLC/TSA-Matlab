% Matlab code for Continous Classification
% Get variable and path inputs
v_input=main_Inputs;
% 1. number of nrows processed
nrows=v_input.ijdim(1);
% 2. number of pixels procesed per line
ncols=v_input.ijdim(2);
% 3. total number of CPU used
tn_cpu=ntasks; 
% 4. locations of data
dir_l=v_input.l_dir;

% idn_cpu => inputs (CPU identification number) for parallel computing
idn_cpu=task;

% load in RF models 
load modelRF500;

% make TSFitMap folder for storing coefficients
n_result=v_input.name_rst;% 'TSFitMap';
if isempty(dir(n_result))
    mkdir(n_result);
end

% Get into the results folder
cd(n_result);

% make TSFitMap folder for storing coefficients
n_record=v_input.name_rec;% 'TSFitMap';
if isempty(dir(n_record))
    mkdir(n_record);
end 

% Get out of the subfolders first
cd ..

% prepare the irows for idn_cpu for ALL rows
irows=zeros(1,1);
i=0;
while idn_cpu+tn_cpu*i<=nrows
   irows(i+1)=idn_cpu+tn_cpu*i;
   i=i+1;
end

for i=1:length(irows)
    % Check whether record_change already exist for row == i
    if isempty(dir([dir_l,'/',n_result,'/','record_change',num2str(irows(i)),'.mat']))
        %          fprintf('Processing the %dth row\n',irows(i));
        
        % Continous Change Detection Done for a line of timeseries pixels
        TrendSeasonalFit_v12Line(dir_l,n_result,n_record,ncols,irows(i),v_input.mini_rmse,...
            v_input.T_cg,v_input.conse,v_input.num_c,v_input.nbands);
    end
    
    % Continous Classfication Done for a Line of timeseries pixels
    Classification_CFT_map(dir_l,n_result,irows(i),modelRF,v_input.num_c,v_input.nbands);
end
exit;