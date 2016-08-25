function ShowAccuChangeMap
% This function is used to provde Change Maps for each pixels
% Revisions: $ Date: 10/19/2015 $ Copyright: Zhe Zhu
% Version 1.1: Add NDVI and NBR for training (10/19/2015)
% Version 1.0: Show change maps (01/10/2015)

addpath('~/ccdc');
v_input = ccdc_Inputs;
v_lct = pwd;

% dimension and projection of the image
nrows = v_input.ijdim(1);
ncols = v_input.ijdim(2);
jiDim = [ncols,nrows];
jiUL = v_input.jiul;
res = v_input.resolu;
zc = v_input.zc;
nbands = v_input.nbands;

% produce confirmed change number map
ChangeNum = zeros(nrows,ncols,1,'uint8'); 
% produce most recent confirmed change time map
ChangeTime=zeros(nrows,ncols,1,'single'); 
% produce most recent confirmed change magnitude map
ChangeMag=zeros(nrows,ncols,1,'single'); 

% % updated change probability map
% PChangeAft = zeros(nrows,ncols,1,'uint8'); 
% % updated probable change time map
% PChangeTimeAft = zeros(nrows,ncols,1,'single');
% % updated probable change magnitude map
% PChangeMagAft = zeros(nrows,ncols,1,'single');
% 
% % updated change probability map
% PChangeBef = zeros(nrows,ncols,1,'uint8'); 
% % updated probable change time map
% PChangeTimeBef = zeros(nrows,ncols,1,'single');
% % updated probable change magnitude map
% PChangeMagBef = zeros(nrows,ncols,1,'single');

% make Predict folder for storing predict images
n_map=v_input.name_map;% 'CCDCMap';
if isempty(dir(n_map))
    mkdir(n_map);
end

% cd to the folder for storing recored structure
cd(v_input.name_rst)

imf=dir('record_change*'); % folder names
num_line=size(imf,1);
line_pt = 0;

for line = 1:num_line
    
    if 100*(line/num_line) - line_pt > 1
        fprintf('Processing %.0f percent\n',ceil(100*(line/num_line)));
        fprintf('%s\n',v_lct);
        line_pt = 100*(line/num_line);
    end
    
    % load one line of time series models
    load(imf(line).name);
    
    % postions 
    pos = [rec_cg.pos];    

    % continue if there is no model available
    l_pos = length(pos);
    if l_pos == 0
        continue;
    end
    
    time = [rec_cg.t_break];
    change_prob=[rec_cg.change_prob];
    magnitude = [rec_cg.magnitude];
    shap_mag = reshape(magnitude,nbands-1,[]);
    % category
    categ = [rec_cg.category];
    
    % for confirmed change
    % ids of confirmed change 
    ids_break = change_prob == 1;
    
    % positions of confirmed change
    pos_break = pos(ids_break);
    % time of confirmed change
    time_break = time(ids_break);
    % magnitude of confirmed change;
    mag_break = shap_mag(:,ids_break);
%     
%     % for probable change
%     % ids of probable change after
%     ids_prob_aft = change_prob > 0 & change_prob < 1;
%     % ids of probable change before
%     ids_prob_bef = change_prob < 0;
%     
%     % positions of probable change
%     pos_prob_aft = pos(ids_prob_aft);
%     pos_prob_bef = pos(ids_prob_bef);
%     % percent probability of probable change
%     pct_prob_aft = change_prob(ids_prob_aft);
%     pct_prob_bef = -change_prob(ids_prob_bef);
%     % time of probable change
%     time_prob_aft = time(ids_prob_aft);
%     time_prob_bef = time(ids_prob_bef);
%     % magnitude of probable change
%     mag_prob_aft = shap_mag(:,ids_prob_aft);
%     mag_prob_bef = shap_mag(:,ids_prob_bef);
    
    % pass if there is no confirmed change
    l_pos = length(pos_break);
    if l_pos > 0  
        for i=1:l_pos
            [I,J] = ind2sub(jiDim,pos_break(i));
            ChangeNum(J,I) = ChangeNum(J,I) + 1;
            vec_time = datevecmx(time_break(i));
            ChangeTime(J,I) = vec_time(1);
            ChangeMag(J,I) = norm(mag_break(:,i));
        end
    end
    
%     % pass if there is no probable change after
%     l_pos = length(pos_prob_aft);
%     if l_pos > 0
%         for i=1:l_pos
%             [I,J] = ind2sub(jiDim,pos_prob_aft(i));
%             PChangeAft(J,I) = uint8(100*pct_prob_aft(i));
%             vec_time = datevecmx(time_prob_aft(i));
%             PChangeTimeAft(J,I) = vec_time(1)*1000 + time_prob_aft(i) - datenummx(vec_time(1),1,0);
%             PChangeMagAft(J,I) = norm(mag_prob_aft(:,i));
%         end
%     end
%     
%     % pass if there is no probable change before
%     l_pos = length(pos_prob_bef);
%     if l_pos > 0
%         for i=1:l_pos
%             [I,J] = ind2sub(jiDim,pos_prob_bef(i));
%             PChangeBef(J,I) = uint8(100*pct_prob_bef(i));
%             vec_time = datevecmx(time_prob_bef(i));
%             PChangeTimeBef(J,I) = vec_time(1)*1000 + time_prob_bef(i) - datenummx(vec_time(1),1,0);
%             PChangeMagBef(J,I) = norm(mag_prob_bef(:,i));
%         end
%     end
end

cd ..
% write ENVI files
% confirmed change
ARD_enviwrite([v_input.l_dir,'/',n_map,'/ChangeNumber'],ChangeNum,'uint8','bsq','example_img');
ARD_enviwrite([v_input.l_dir,'/',n_map,'/ChangeYear'],ChangeTime,'single','bsq','example_img');
ARD_enviwrite([v_input.l_dir,'/',n_map,'/ChangeMagnitude'],ChangeMag,'single','bsq','example_img');
% 
% % probable change after
% enviwrite([v_input.l_dir,'/',n_map,'/PChangeProbAft'],PChangeAft,'uint8',res,jiUL,'bsq',zc);
% enviwrite([v_input.l_dir,'/',n_map,'/PChangeTimeAft'],PChangeTimeAft,'single',res,jiUL,'bsq',zc);
% enviwrite([v_input.l_dir,'/',n_map,'/PChangeMagAft'],PChangeMagAft,'single',res,jiUL,'bsq',zc);
% 
% % probable change before
% enviwrite([v_input.l_dir,'/',n_map,'/PChangeProbBef'],PChangeBef,'uint8',res,jiUL,'bsq',zc);
% enviwrite([v_input.l_dir,'/',n_map,'/PChangeTimeBef'],PChangeTimeBef,'single',res,jiUL,'bsq',zc);
% enviwrite([v_input.l_dir,'/',n_map,'/PChangeMagBef'],PChangeMagBef,'single',res,jiUL,'bsq',zc);