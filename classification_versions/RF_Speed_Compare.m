% function overall = RF_Speed_Compare(N)
% Function for accuarcy assessment
% Design different strategies (Equal number)

% CCDC 1.6 version - Zhe Zhu, EROS, USGS
%
% Revisions: $ Date: 05/11/2015 $ Copyright: Zhe Zhu
% Version 1.6  Change back to original classification method (05/11/2015)
% Version 1.5  Only use undisturbed data (05/11/2015)
% Version 1.4  Hide one block and do accuracy assessment (03/02/2015)
% Version 1.3: Use version 7.3 for storing RF model (01/10/2015)
% Version 1.2: Add ancillary data from NLCD (01/10/2015)
% Version 1.1: Fixed a bug in picking the wrong pixel for training (11/08/2014)

% Inputs:
% image with DN values show what land cover class (Trends)
% image with DN values show what block the data are from (Trends_ids)

% Use the default fmask toobox developed by Zhe Zhu
addpath('~/ccdc');
% Tools of RFC
addpath('~/Algorithms/CCDC/Tools/RF_Class_C');

%% Equal proportion 
clear
clc
pwd
% get input data
load('Ys_up');
load('Xs_up');

% repeat number
N = 1;
n_rep = N;

% remove Y == 3 or 10 disturbed classes
ids_rm = Y(:,1) == 3 | Y(:,1) == 10 | Y(:,1) == 0 | Y(:,2) == 0;
X(ids_rm,:) = [];
Y(ids_rm,:) = [];
% Tradditional inputs
X = X(:,1:7*9);

% get class number
all_class = unique(Y(:,1));
% number of class
n_class = length(all_class);
% equal rate (0.01%, 0.02% ... 0.1%)
eq_num = 10000; % optimum #
% eq_num = [18000 20000];
% number of test scenarios
n_eq = length(eq_num);
% k-folder validation
k_v = unique(Y(:,2))';
% % accuracies
% rf_accuracy = zeros(n_eq,n_rep);
% % weighted accuracies
% wt_accuracy = zeros(n_eq,n_rep);
% % minimum accuracies
% mi_accuracy = zeros(n_eq,n_rep);
% result matrix
ccdc_res = zeros(n_eq,n_rep,4);

% calculate proportion based # for each class
prct = hist(Y(:,1),all_class);
prct = prct/sum(prct);

overall = zeros(N*length(k_v),3);

for i_rep = 1:n_rep
    for i_eq = 1:n_eq
        % confusion matrix
        % confm = zeros(n_class,n_class);
        % y vs yhat
        Y_hat_c = [];
        Y_tst_c = [];
        
        % k-folder cross validation
        for i_q = k_v
            % get ids of validation and training data
            id_trn = Y(:,2) ~= i_q;
            id_tst = Y(:,2) == i_q;
            
            % get training data based on basic CCDC inputs
            X_trn = X(id_trn,:);
            Y_trn = Y(id_trn,1);
            
            % testing data based on basic CCDC inputs
            X_tst = X(id_tst,:);
            Y_tst = Y(id_tst,1);
            
            % intialized selected X & Y training data
            sel_X_trn = [];
            sel_Y_trn = [];
            
            for i_class = 1:n_class
                % find ids for each land cover class
                ids = find(Y_trn == all_class(i_class));
                % total # of reference pixels for permute
                tmp_N = length(ids);
                
                % random permute the reference pixels
                tmp_rv = randperm(tmp_N);
                
                % adjust num_prop based on proportion
                adj_num = ceil(eq_num(i_eq)*prct(i_class));
                
                if tmp_N > adj_num
                    % if tmp_N > adj_num, use adj_num, otherwise, use tmp_N
                    tot_n = adj_num;
                else
                    tot_n = tmp_N;
                end
                
                % permutted ids
                rnd_ids = ids(tmp_rv(1:tot_n));
                
                % X_trn and Y_trn
                sel_X_trn = [sel_X_trn; X_trn(rnd_ids,:)];
                sel_Y_trn = [sel_Y_trn; Y_trn(rnd_ids)];
            end
            
            tic
            for i=1:N
                model = classRF_train(sel_X_trn,sel_Y_trn,100);
                Y_hat = classRF_predict(X_tst,model);
                overall((i_q-1)*N+i:i_q*N,1) = 100*sum(Y_hat==Y_tst)/length(Y_tst)
            end
            toc
            
            tic
            for i=1:N
                B = fitensemble(sel_X_trn,sel_Y_trn,'Bag',100,'Tree','type','classification');
                Y_hat = (B.predict(X_tst));
                overall((i_q-1)*N+i:i_q*N,2) = 100*sum(Y_hat==Y_tst)/length(Y_tst)
            end
            toc
            
            tic
            for i=1:N
                B = fitensemble(sel_X_trn,sel_Y_trn,'AdaBoostM2',100,'Tree','type','classification');
                Y_hat = (B.predict(X_tst));
                overall((i_q-1)*N+i:i_q*N,3) = 100*sum(Y_hat==Y_tst)/length(Y_tst)
            end
            toc                       
        end
        
    end
end
fprintf('Fishined Equal Number!\n');