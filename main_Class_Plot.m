function [map,qa,rec_cg] = main_Class_Plot(row,col,yy,mm,dd)
%
% Revisions: $ Date: 07/03/2015 $ Copyright: Zhe Zhu
% Version 1.2: Classify based on best strategy (07/03/2015)
% Version 1.1: Add ancillary data from NLCD & Fmask (01/10/2015)
%
% addpath('~/ccdc/TSFit_versions/');
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
load modelRF1_6 % change name if neccessary

% load ancillary data (A)
n_anc = 5 + 3; % NLCD (5) + Fmask (3)
anc = zeros(ncols,nrows,n_anc);
% get to the reference and anc data folder
cd ANC;

im_tmp = double(enviread('aspect')); % 1.dem deriv
anc(:,:,1) = im_tmp';
im_tmp = double(enviread('dem')); % 2. dem deriv
anc(:,:,2) = im_tmp';
im_tmp = double(enviread('posidex')); % 3. dem deriv
anc(:,:,3) = im_tmp';
im_tmp = double(enviread('slope')); % 4. dem deriv
anc(:,:,4) = im_tmp';
im_tmp = double(enviread('mpw')); % 5. wetland potential index
anc(:,:,5) = im_tmp';

im_tmp = double(enviread('Fmask_stat')); % 6~8. Fmask statistics
anc(:,:,6) = im_tmp(:,:,1)';
anc(:,:,7) = im_tmp(:,:,2)';
anc(:,:,8) = im_tmp(:,:,3)';
cd ..

map = zeros(1,length(yy));
qa = zeros(1,length(yy));
% Continous Classfication Done for a Line of timeseries pixels
rec_cg = TrendSeasonalFit_v12Plot(row,col,0.99,6,2:6);

for i = 1:length(yy)
    [map(i),qa(i)] = Classification_Plot(rec_cg,anc,modelRF,yy(i),mm,dd);
end

map = uint16([map;yy]);
qa = uint16([qa;yy]);
% class = map;
% % whether it is a transition class
% if class(1) == class(end)
%     % no transition
%     class(:) = class(end);
% else
%     % trainsition class
%     id_change = find(class~=class(1));
%     % update classes based on transition id
%     class(1:id_change(1)-1) = class(1);
%     class(id_change(1):end) = class(end);
% end
end
