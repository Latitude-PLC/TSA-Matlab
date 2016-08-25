function [map_class,map_pos] = Classification_CFT_map_Temp(dir_l,n_rst,nrow,model,num_c,nbands,anc,map_t)
% This is the classification algorithm for classifying all map pixels by lines

% Revisions: $ Date: 02/10/2015 $ Copyright: Zhe Zhu
% Version 1.2: Add in the temporal information to the training (02/10/2015)
% Version 1.1: Add ancillary data from NLCD (01/10/2015)

% initiallize output
map_class = [];
map_pos = [];

% fprintf('Processing the %d row\n',i);
load([dir_l,'/',n_rst,'/','record_change',num2str(nrow)]);

% get information for classification
pos=[rec_cg.pos]; % position
l_pos = length(pos);

if l_pos > 0
    t_start = [rec_cg.t_start];% start time
    rmse=[rec_cg.rmse]; % rmse
    coefs = [rec_cg.coefs]; % coefficients
    % reshape coefs
    coefs = reshape(coefs,num_c,nbands-1,[]);
    
    % curve id that are exist before map_tget
    id_lgt = map_t > t_start;
    pos = pos(id_lgt);
    t_start = t_start(id_lgt);
    rmse = rmse(:,id_lgt);
    coefs = coefs(:,:,id_lgt);
    
    % curve ids that the last for the same location
    [pos,id_last] = unique(pos,'last');
    t_start = t_start(id_last);
    rmse = rmse(:,id_last);
    coefs = coefs(:,:,id_last);
    
    % number of curves per line
    l_pos=length(pos);
    % number of bands for ancillary data (A)
    n_anc = size(anc,3);
    
    if l_pos > 0 % has more than one curves exsit for each line
        % Prepare for classification inputs
        Xclass=zeros(l_pos,(num_c+1)*(nbands-1)+n_anc);
        
        % get ancillary data (B)
        % array for storing ancillary data
        array_anc = zeros(n_anc,l_pos);
        for i_a = 1:n_anc - 1
            tmp_anc = anc(:,:,i_a);
            array_anc(i_a,:) = tmp_anc(pos);
        end
        
        for icol = 1:l_pos
            % modified constant as inputs
            coefs(1,:,icol) = coefs(1,:,icol)+map_t*coefs(2,:,icol);
            cft_tmp = coefs(:,:,icol);
            % add rmse & anc for each band as  (C)
            % add model time to the last band of anc
            array_anc(end,icol) = map_t - t_start(icol);
            Xclass(icol,:)=[rmse(:,icol);cft_tmp(:);array_anc(:,icol)];
        end
        
        map_class = classRF_predict(Xclass,model); % class
        map_pos = pos; % position
        
    end % end of if num_s > 0
end % end of if l_pos > 0

end % end of the function

