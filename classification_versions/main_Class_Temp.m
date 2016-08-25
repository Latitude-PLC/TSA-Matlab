% Matlab code for Continous Classification
%
% Revisions: $ Date: 01/10/2015 $ Copyright: Zhe Zhu
% Version 1.1: Add ancillary data from NLCD (01/10/2015)
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
load('/projectnb/landsat/projects/LCMAP/p046r027/images/v12.17.CCDC/v12.17.Model/modelRF_50k_1k_Temp.mat'); % change name if neccessary
%% load ancillary data (A)
% n_anc = 6 + 1 + 4 + 1; % nlcd + cgnum + fmask + changed time
% anc = zeros(ncols,nrows,n_anc);
% % get to the reference and anc data folder
% dir_anc = '/ccdc_train/anc/';
% im_tmp = double(enviread([dir_l,dir_anc,'aspect'])); % 1.dem deriv
% anc(:,:,1) = im_tmp';
% im_tmp = double(enviread([dir_l,dir_anc,'cti'])); % 2. dem deriv
% anc(:,:,2) = im_tmp';
% im_tmp = double(enviread([dir_l,dir_anc,'dem'])); % 3. dem deriv
% anc(:,:,3) = im_tmp';
% im_tmp = double(enviread([dir_l,dir_anc,'posidex'])); % 4. dem deriv
% anc(:,:,4) = im_tmp';
% im_tmp = double(enviread([dir_l,dir_anc,'slope'])); % 5. dem deriv
% anc(:,:,5) = im_tmp';
% im_tmp = double(enviread([dir_l,dir_anc,'wpi'])); % 6. wetland potential index
% anc(:,:,6) = im_tmp';
% 
% 
% im_tmp = double(enviread([dir_l,dir_anc,'ChangeNum'])); % 7. change #
% anc(:,:,7) = im_tmp';
% 
% im_tmp = double(enviread([dir_l,dir_anc,'Fmask_stat'])); % 8-11. Fmask statistics
% anc(:,:,8) = im_tmp(:,:,1)';
% anc(:,:,9) = im_tmp(:,:,2)';
% anc(:,:,10) = im_tmp(:,:,3)';
% anc(:,:,11) = im_tmp(:,:,4)';
% % last band fill with zeros

% load ancillary data (A)
% n_anc = 6 + 1 + 4 + 1; % nlcd + cgnum + fmask + changed time
n_anc = 1 + 3; % dem + fmask (1,3,4)
anc = zeros(ncols,nrows,n_anc);
% get to the reference and anc data folder
dir_anc = '/ANC/';
im_tmp = double(enviread([dir_l,dir_anc,'dem'])); % 3. dem deriv
anc(:,:,1) = im_tmp';
im_tmp = double(enviread([dir_l,dir_anc,'Fmask_stat'])); % 8-11. Fmask statistics
anc(:,:,2) = im_tmp(:,:,1)'; % land 
anc(:,:,3) = im_tmp(:,:,3)'; % snow
anc(:,:,4) = im_tmp(:,:,4)'; % cloud

% make TSFitMap folder for storing coefficients
n_result=v_input.name_rst;% 'TSFitMap';
if isempty(dir(n_result))
    mkdir(n_result);
end

% prepare the irows for task for ALL rows
irows=zeros(1,1);
i=0;
while task+ntasks*i<=nrows
   irows(i+1)=task+ntasks*i;
   i=i+1;
end

for i = 1:length(irows)
    
    % Check whether record_change already exist for row == i
    if isempty(dir([dir_l,'/',n_result,'/','record_change',num2str(irows(i)),'.mat']))                
        % Continous Change Detection Done for a line of timeseries pixels
        TrendSeasonalFit_v12Line(dir_l,n_result,n_record,ncols,irows(i),v_input.mini_rmse,...
            v_input.T_cg,v_input.conse,v_input.num_c,v_input.nbands);
    end
    
    % background name for identification of row (0~9999)
    id_name='0000';
    % str number for identification of row
    str_num=num2str(irows(i));
    % add str number to background
    id_name((end-length(str_num)+1):end) = str_num;
    
%    if isempty(dir([dir_l,'/',n_result,'/','TSFitMapMat',id_name,'.mat']))
        fprintf('Processing the %dth row\n',irows(i));
        % Continous Classfication Done for a Line of timeseries pixels (B)
        Classification_CFT_map(dir_l,n_result,irows(i),modelRF,v_input.num_c,v_input.nbands,anc);
%    end
end
fprintf('Done!\n');
exit;