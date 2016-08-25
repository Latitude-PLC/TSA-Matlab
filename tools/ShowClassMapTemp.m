% function ShowClassMapTemp(yy,mm,dd)
% This function is used to provde Change Type Maps for each pixels
yy=2000;
mm = 6;
dd = 1;
[yy mm dd] %#ok<NOPTS>
addpath('~/ccdc');
v_input = ccdc_Inputs;
pwd
% nbands=v_input.nbands-1;% number of bands in the image (except the Fmask)
% ncoefs=v_input.num_c;% number of coefficients

% dimension and projection of the image
nrows = v_input.ijdim(1);
ncols = v_input.ijdim(2);
jiDim = [ncols,nrows];
jiUL = v_input.jiul;
res = v_input.resolu;
zc = v_input.zc;
dir_l = v_input.l_dir;
num_c = v_input.num_c; 
nbands = v_input.nbands;

% converted to julian date
j_date=datenum(yy,mm,dd);
% date for show i.e. 19990815
s_date=yy*10000+mm*100+dd;

% load in RF models 
load modelRF_50k_1k_Temp_ANC; % change name if neccessary

% load ancillary data (A)
n_anc = 1 + 3 + 1 + 1; % dem + fmask (1,3,4) + cg# + cgtime
anc = zeros(ncols,nrows,n_anc);
% get to the reference and anc data folder
dir_anc = '/ANC/';
im_tmp = double(enviread([dir_l,dir_anc,'dem'])); % 3. dem deriv
anc(:,:,1) = im_tmp';
im_tmp = double(enviread([dir_l,dir_anc,'Fmask_stat'])); % 8-11. Fmask statistics
anc(:,:,2) = im_tmp(:,:,1)'; % land 
anc(:,:,3) = im_tmp(:,:,3)'; % snow
anc(:,:,4) = im_tmp(:,:,4)'; % cloud
im_tmp = double(enviread([dir_l,dir_anc,'ChangeNum'])); % 7. change #
anc(:,:,5) = im_tmp'; % change #

% produce land cover map at any date
ClassType = zeros(nrows,ncols,1,'uint8'); 
% transpose matrix
tClassType = ClassType';

% make Predict folder for classification maps 
n_map = v_input.name_map;% 'CCDCMap';
if isempty(dir(n_map))
    mkdir(n_map);
end

% make TSFitMap folder for storing coefficients
n_result=v_input.name_rst;% 'TSFitMap';
if isempty(dir(n_result))
    mkdir(n_result);
end

% initialize process percentage
pro_pct = -10;
for i_row = 1:nrows
    % printf status
    if round(100*(i_row/nrows)) - pro_pct >= 10
        pro_pct = round(100*(i_row/nrows));
        fprintf('%d%%...',pro_pct);
    end 

    [class,pos] = Classification_CFT_map_Temp(dir_l,n_result,i_row,modelRF,num_c,nbands,anc,j_date);
    % label class 
    tClassType(pos) = class;
    
    % end of process
    if i_row == nrows
        fprintf('Done.\n');
    end
end

% tranpose back
ClassType = tClassType';

% write ENVI files
enviwrite([dir_l,'/',n_map,'/Trends_Temporal_ANC',...
    num2str(s_date)],ClassType,'uint8',res,jiUL,'bsq',zc);