% Function for extracting training data
% Prepare data for Training_Strategy.m
%
% CCDC 1.3 version - Zhe Zhu, EROS, USGS
%
% Revisions: $ Date: 11/25/2015 $ Copyright: Zhe Zhu
% Version 1.3  Select sample based on best strategy (06/30/2015)
% Version 1.2  Only use undisturbed data (05/11/2015)
% Version 1.1  Use version 7.3 for storing RF model (01/10/2015)
% Version 1.0  Fixed a bug in picking the wrong pixel for training (11/08/2014)
%
% Inputs:
% image with DN values show what land cover class (Trends)
% image with DN values show what block the data are from (Trends_ids)

%% Prepare for the inputs
clear
clc
v_input = ccdc_Inputs;
pwd
% Inputs:
% 1. Location of "im_roi" & place to save "modelRF"
l_dir = v_input.l_dir;

% 2. ground truth time interval
gt_start = v_input.gt(1);
gt_end = v_input.gt(2);

% Constants:
% number of coefficients
num_c = v_input.num_c;
% number of bands
nbands = v_input.nbands; 
% image dimension
nrows = v_input.ijdim(1);
ncols = v_input.ijdim(2);
jiDim = [ncols,nrows];
% version of CCDC
ccdc_v = 1.3;

% use the Trends
im_roi = enviread('example_img'); % Land Cover Trends data (all pixels)

% transformed to Matlab dimension (opposite in x & y)
tim_roi = im_roi'; 
clear im_roi;

% Find nonzero ids with 1 2 3
%                       4 5 6 
idsfind = find(tim_roi > 0);
[~,i_ids] = ind2sub(jiDim,idsfind);

% length of roi pixels
rec_l = length(idsfind);

% Get Training data prepared
% the training data is NxD and labels are Nx1, where N=#of
% examples, D=# of features
% intialize maximum Xs and Ys
X = zeros(rec_l,(num_c+1)*(nbands-1)); % 7 bands cft & rmse
Y = zeros(rec_l,2); % Trends classes + location

% Get into the TSFitMap folder
cd(v_input.name_rst);
% intiate i_row
i_row = -1;
% number of pixels for traning
plusid = 0;

for i=1:rec_l
    % Just load once for a line of rec_cg for all reference within this line    
    if i_ids(i) ~= i_row
        fprintf('Processing the %dth line ...\n',i_row);
        % load CCDCRec
        load(['record_change',num2str(i_ids(i))]);
        
        % matrix of each component
        t_start = [rec_cg.t_start];
        t_end = [rec_cg.t_end];
        coefs = [rec_cg.coefs];
        rmse = [rec_cg.rmse];
        pos = [rec_cg.pos];
        categ = [rec_cg.category];
        % reshape coefs
        coefs = reshape(coefs,num_c,nbands-1,[]);       
    end
       
    % find the curve within a fixed time interval
    ids_line = find(pos == idsfind(i));
    
    for j = 1:length(ids_line)
        % id of reference data
        id_ref = ids_line(j);
        % position of reference data
        pos_ref = pos(id_ref);
        
        % take curves that fall witin the training period
        % remove curves that are changed within training period
        if t_start(id_ref) < gt_start & t_end(id_ref) > gt_end 
            % number of time series model for training
            plusid = plusid + 1;
            
            % two way of overall reflectance
            tmp_cft = coefs(:,:,id_ref);
            % temporal overall ref
            % tmp_cft(1,:) = tmp_cft(1,:)+gt_mid*tmp_cft(2,:);
            tmp_cft(1,:) = tmp_cft(1,:) + 0.5*(t_start(id_ref)+t_end(id_ref))*tmp_cft(2,:);
          
            % all rmse
            tmp_rmse = rmse(:,id_ref);            

            % prepare Xs
            X(plusid,:) = [tmp_rmse;tmp_cft(:)];
            Y(plusid,1) = tim_roi(pos_ref); % class category
            Y(plusid,2) = pos_ref;
            % Y(plusid,2) = tim_bw(pos_ref); % block #
            % Y(plusid,3) = tim_roi(pos_ref); % block #
        end
    end
    i_row = i_ids(i);        
end

% remove out of boundary or changed pixels
X = X(1:plusid,:);
Y = Y(1:plusid,:);

% make sure they are double!
Y = double(Y);
X = double(X);

% cd to the images folder
cd(l_dir);

% save Xs and Ys
save('Ys','Y','-v7.3');
save('Xs','X','-v7.3');

%% selecting pixels for training
% get input data
load('Ys');
load('Xs');

% only use the first column for training
Y = Y(:,1);

% remove Y == 3 or 10 disturbed classes
ids_rm = Y == 3 | Y == 10 | Y == 0;
X(ids_rm,:) = [];
Y(ids_rm) = [];

% number of variables
x_dim = size(X,2);

% update class number
all_class = unique(Y);
% update number of class
n_class = length(all_class);
% calculate proportion based # for each class
prct = hist(Y(:,1),all_class);
prct = prct/sum(prct);

% number of reference for euqal number training
eq_num = 20000; % total # 
n_min = 600; % minimum # 
n_max = 8000; % maximum # 

% intialized selected X & Y training data
sel_X_trn = [];
sel_Y_trn = [];

for i_class = 1:n_class
    % find ids for each land cover class
    ids = find(Y == all_class(i_class));
    % total # of reference pixels for permute
    tmp_N = length(ids);
    
    % random permute the reference pixels
    tmp_rv = randperm(tmp_N);
    
    % adjust num_prop based on proportion
    adj_num = ceil(eq_num*prct(i_class));
    
    % adjust num_prop based on min and max
    if adj_num < n_min
        adj_num = n_min;
    elseif adj_num > n_max
        adj_num = n_max;
    end
    
    if tmp_N > adj_num
        % if tmp_N > adj_num, use adj_num, otherwise, use tmp_N
        tot_n = adj_num;
    else
        tot_n = tmp_N;
    end
    
    % permutted ids
    rnd_ids = ids(tmp_rv(1:tot_n));
    
    % X_trn and Y_trn
    sel_X_trn = [sel_X_trn; X(rnd_ids,:)];
    sel_Y_trn = [sel_Y_trn; Y(rnd_ids)];
end

% log for CCDC Train paramters and versions
% report only for the first task
class_value = unique(sel_Y_trn);
class_number = hist(sel_Y_trn,class_value);

fileID = fopen('CCDC_Train_log.txt','w');
% write location of image stack
fprintf(fileID,'Image location = %s\r\n',l_dir);
% write number of images used
fprintf(fileID,'Number of sample = %d\r\n',sum(class_number));
% write number of images used
fprintf(fileID,'Minimum number of sample per class = %d\r\n',min(class_number));
% write number of images used
fprintf(fileID,'Maximum number of sample per class = %d\r\n',max(class_number));
% CCDC Version
fprintf(fileID,'CCDC Train based on Trends Version = %.2f\r\n',ccdc_v);
% updates
fprintf(fileID,'******************************************************************************************************\r\n');
fprintf(fileID,'Revisions: $ Date: 11/20/2015 $ Copyright: Zhe Zhu\r\n');
fprintf(fileID,'Version 1.3   Select sample based on best strategy (06/30/2015)\r\n');
fprintf(fileID,'Version 1.2   Only use undisturbed data (05/11/2015)\r\n');
fprintf(fileID,'Version 1.1   Use version 7.3 for storing RF model (01/10/2015) \r\n');
fprintf(fileID,'Version 1.0   Fixed a bug in picking the wrong pixel for training (11/08/2014)\r\n');
fprintf(fileID,'******************************************************************************************************\r\n');
fclose(fileID);

modelRF = classRF_train(sel_X_trn,sel_Y_trn);
save('modelRF','modelRF','-v7.3'); 
fprintf('Finishined Training!\n');

