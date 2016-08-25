function Map = Classification_CFT_map(dir_l,n_rst,nrow,model,num_c,nbands,anc)
% This is the classification algorithm for classifying all map pixels by lines
%
% Revisions: $ Date: 07/03/2015 $ Copyright: Zhe Zhu
% Version 1.2: Classify based on best strategy (07/03/2015)
% Version 1.1: Add ancillary data from NLCD (01/10/2015)

% fprintf('Processing the %d row\n',i);
load([dir_l,'/',n_rst,'/','record_change',num2str(nrow)]);
% load data
% Map1 (t_start) 
% Map2 (t_end) 
% Map3 (pos) 
% Map4 (class) 
% Map5 (t_break) 
% rec_cgMap6 (change_prob)

t_start=[rec_cg.t_start];% Map 1
t_end=[rec_cg.t_end]; % Map 2
pos=[rec_cg.pos]; % Map 3
% class type = Map 4
change_prob = [rec_cg.change_prob]; % Map 5
t_break=[rec_cg.t_break];% Map 6
rmse=[rec_cg.rmse]; % each curve has nbands-1 rmse
% number of curves per line
num_s=length(pos);
%% number of bands for ancillary data (A)
n_anc = size(anc,3);
%% end of (A)

if num_s > 0 % has more than one curves exsit for each line
    % model coefficients
    tmp=[rec_cg.coefs];
    % Prepare for classification inputs
    Xclass=zeros(num_s,(num_c+1)*(nbands-1)+n_anc);
    
    % initialaze Map
    Map=zeros(5,num_s);
    Map(1,:)=t_start;
    Map(2,:)=t_end;
    Map(3,:)=pos;
    
    %% get ancillary data (B)
    % array for storing ancillary data
    array_anc = zeros(n_anc,num_s);
    for i_a = 1:n_anc%-1
        tmp_anc = anc(:,:,i_a);
        array_anc(i_a,:) = tmp_anc(pos);
    end   
    %% end of (B)
    
    for icol=1:num_s;
        % coefficients from the 7 bands
        i_tmp=tmp(:,((icol-1)*(nbands-1)+1):(nbands-1)*icol);
        % modified constant as inputs
        i_tmp(1,:)=i_tmp(1,:)+(t_start(icol)+t_end(icol))*i_tmp(2,:)/2;
        %% add rmse & anc for each band as  (C)
        % % add model time to the last band of anc
        % array_anc(end,icol) = t_start(icol) - t_end(icol);
        Xclass(icol,:)=[reshape(rmse(((icol-1)*(nbands-1)+1):(nbands-1)*icol),nbands-1,1);i_tmp(:);array_anc(:,icol)];
        %% end of (C)
    end
    
    Map(4,:)=classRF_predict(Xclass,model); % class
    Map(5,:)=change_prob; % change probability
    Map(6,:)=t_break; % time of breark
    
    % write TSFitMapMat.mat
    % background name for identification of row (0~9999)
    id_name='0000';
    % str number for identification of row
    str_num=num2str(nrow);
    % add str number to background
    id_name((end-length(str_num)+1):end)=str_num;
    save([dir_l,'/',n_rst,'/','TSFitMapMat',id_name],'Map');
end % end of if

end % end of the function

