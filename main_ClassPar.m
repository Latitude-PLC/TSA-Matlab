% Matlab code for Continous Classification
%
% Revisions: $ Date: 10/19/2015 $ Copyright: Zhe Zhu
% Version 1.3: Add NDVI and NBR for training (10/19/2015)
% Version 1.2: Classify based on best strategy (07/03/2015)
% Version 1.1: Add ancillary data from NLCD & Fmask (01/10/2015)
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

% load in RF models 
load modelRF % change name if neccessary

% load ancillary data (A)
n_anc = 5 + 3; % NLCD (5) + Fmask (3)
anc = zeros(ncols,nrows,n_anc);
% get to the reference and anc data folder
dir_anc = '/ANC/';

im_tmp = double(enviread([dir_l,dir_anc,'aspect'])); % 1.dem deriv
anc(:,:,1) = im_tmp';
im_tmp = double(enviread([dir_l,dir_anc,'dem'])); % 2. dem deriv
anc(:,:,2) = im_tmp';
im_tmp = double(enviread([dir_l,dir_anc,'posidex'])); % 3. dem deriv
anc(:,:,3) = im_tmp';
im_tmp = double(enviread([dir_l,dir_anc,'slope'])); % 4. dem deriv
anc(:,:,4) = im_tmp';
im_tmp = double(enviread([dir_l,dir_anc,'mpw'])); % 5. wetland potential index
anc(:,:,5) = im_tmp';

im_tmp = double(enviread([dir_l,dir_anc,'Fmask_stat'])); % 6~8. Fmask statistics
anc(:,:,6) = im_tmp(:,:,2)';
anc(:,:,7) = im_tmp(:,:,3)';
anc(:,:,8) = im_tmp(:,:,4)';

% make TSFitMap folder for storing coefficients
n_result=v_input.name_rst;% 'TSFitMap';
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
    
%     % Check whether record_change already exist for row == i
%     if isempty(dir([dir_l,'/',n_result,'/','record_change',num2str(irows(i)),'.mat']))
%         fprintf('Missing the %dth row!\n',irows(i));
%         continue
%     end
%     
%     fprintf('Processing the %dth row\n',irows(i));
%     % Continous Classfication Done for a Line of timeseries pixels
%     Classification_Line(dir_l,n_result,irows(i),modelRF,v_input.num_c,v_input.nbands,anc);

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
    Classification_Line(dir_l,n_result,irows(i),modelRF,v_input.num_c,v_input.nbands,anc,500);
end
fprintf('Done!\n');
exit;