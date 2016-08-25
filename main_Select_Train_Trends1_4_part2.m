% Function for training data
% Prepare data for Training_Strategy.m
%
% CCDC 1.4 version - Zhe Zhu, EROS, USGS
%
% Revisions: $ Date: 11/25/2015 $ Copyright: Zhe Zhu
%
% Version 1.4  Add ancillary data from NLCD and Fmask (11/25/2015)
% Version 1.3  Select sample based on best strategy (06/30/2015)
% Version 1.2  Only use undisturbed data (05/11/2015)
% Version 1.1  Use version 7.3 for storing RF model (01/10/2015)
% Version 1.0  Fixed a bug in picking the wrong pixel for training (11/08/2014)
%
% Inputs:
% image with DN values show what land cover class (Trends)
% image with DN values show what block the data are from (Trends_ids)

%% selecting pixels for training
function main_Select_Train_Trends1_4_part2(n_times)

% get input data
load('Ys');
load('Xs');

% number of times the area of a standard Landsat scene
n_times = n_times*(25/37);

% only use the first column for training
Y = Y(:,1);  %#ok<NODEF>

% combine mining (4) into developed (2)
Y(Y == 4) = 2;

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
eq_num = ceil(20000*n_times); % total # 
n_min = ceil(600*n_times); % minimum # 
n_max = ceil(8000*n_times); % maximum # 

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

modelRF = classRF_train(sel_X_trn,sel_Y_trn);
save('modelRF','modelRF','-v7.3'); 
fprintf('Finishined Training!\n');
end

