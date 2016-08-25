function [class,classQA] = Classification_Plot(rec_cg,anc,model,yy,mm,dd)
% This is the classification algorithm for classifying all map pixels by lines
%
% Revisions: $ Date: 07/05/2015 $ Copyright: Zhe Zhu
% Version 1.3: Combine classification and change detection output (07/05/2015)
% Version 1.2: Classify based on best strategy (07/03/2015)
% Version 1.1: Add ancillary data from NLCD (01/10/2015)

% pwd

num_c = 8;
nbands = 8;
class = 0;
classQA = 0;

% total number of time series models (including empty ones)
num_s =  length(rec_cg); 
% number of trees
ntrees = 500;

% start time
t_start = [rec_cg.t_start];
% end time
t_end = [rec_cg.t_end];
% position
pos = [rec_cg.pos];
% rmse
rmse = [rec_cg.rmse]; % each curve has nbands-1 rmse
% number of bands for ancillary data
n_anc = size(anc,3);

% model coefficients
tmp = [rec_cg.coefs];
% prepare for classification inputs
Xclass = zeros(1,(num_c + 1)*(nbands - 1) + n_anc);

% find the curve for the specific time intervals
t_cv = datenum(yy,mm,dd);
icol = find(t_start < t_cv & t_end > t_cv);

% array for storing ancillary data
array_anc = anc(pos(2),pos(1),:);

if sum(icol) > 0   
    % coefficients from the 7 bands
    i_tmp = tmp(:,((icol-1)*(nbands-1)+1):(nbands-1)*icol);
    i_rmse = reshape(rmse(((icol-1)*(nbands-1)+1):(nbands-1)*icol),nbands-1,1);
%     % test the slope
%     test = i_tmp(2,:)*(t_end(icol)-t_start(icol)+1)./i_rmse';
%     
%     test_v = norm(test(2:6))^2;
%     T_cg = chi2inv(0.99,5);
%     
%     if test_v > T_cg
%         fprintf('test of dynamics = %.4f\n',test_v);
%     else
%         fprintf('stable class = %.4f\n',test_v);
%     end
    
    % modified constant as inputs
    % i_tmp(1,:) = i_tmp(1,:)+(t_start(icol)+t_end(icol))*i_tmp(2,:)/2;
    i_tmp(1,:) = i_tmp(1,:) + t_cv*i_tmp(2,:);
    % input ready!
    Xclass(:) = [i_rmse;i_tmp(:);array_anc(:)];
    
    % classify the whole line
    [class,votes] = classRF_predict(Xclass,model,ntrees); % class    
else
    class = 0;
    votes = 0;
end

% largest number of votes
[max_v1,max_id] = max(votes(:));
% make this smallest
votes(max_id) = 0;
% second largest number of votes
max_v2 = max(votes(:));
% provide unsupervised ensemble margin as QA
classQA = 100*(max_v1-max_v2)/ntrees;


end % end of the function

