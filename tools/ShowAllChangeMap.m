function ShowAllChangeMap
% This function is used to provde all change maps for each year
% Results for LCMS
% INPUTS:
% % maximum number of images show change time and magnitude
max_n = 30; % 1985~2015

v_input=main_Inputs;
pwd

% dimension and projection of the image
nrows = v_input.ijdim(1);
ncols = v_input.ijdim(2);
jiDim = [ncols,nrows];
jiUL = v_input.jiul;
res = v_input.resolu;
zc = v_input.zc;
nbands = v_input.nbands-1;

% produce confirmed change number map
ChangeNum = zeros(nrows,ncols,1,'uint8'); 
% produce most recent confirmed change time map
ChangeTime=zeros(nrows,ncols,max_n,'uint16'); 
% produce most recent confirmed change magnitude map
ChangeMag=zeros(nrows,ncols,max_n,'single'); 

% updated change probability map
PChangeProb = zeros(nrows,ncols,1,'uint8'); 
% updated probable change time map
PChangeTime = zeros(nrows,ncols,1,'uint16');
% updated probable change magnitude map
PChangeMag = zeros(nrows,ncols,1,'single');

% change probability map at the begining
PChangeProbBef = zeros(nrows,ncols,1,'uint8'); 
% updated probable change time map
PChangeTimeBef = zeros(nrows,ncols,1,'uint16');
% updated probable change magnitude map
PChangeMagBef = zeros(nrows,ncols,1,'single');

% make Predict folder for storing predict images
n_map=v_input.name_map;% 'CCDCMap';
if isempty(dir(n_map))
    mkdir(n_map);
end

% cd to the folder for storing recored structure
cd(v_input.name_rst);

imf=dir('record_change*'); % folder names
num_line=size(imf,1);
for line=1:num_line
    fprintf('Processing %.2f percent\n',100*(line/num_line));
    load(imf(line).name);
    
    % postions 
    pos = [rec_cg.pos];

    % continue if there is no model available
    l_pos=length(pos);
    if l_pos==0
        continue;
    end
    time = [rec_cg.t_break];
    change_prob=[rec_cg.change_prob];
    magnitude = [rec_cg.magnitude];
    shap_mag = reshape(magnitude,nbands,[]);
    
    % for confirmed change
    % ids of confirmed change
    ids_break = change_prob == 1;
    
    % positions of confirmed change
    pos_break = pos(ids_break);
    % time of confirmed change
    time_break = time(ids_break);
    % magnitude of confirmed change;
    mag_break = shap_mag(:,ids_break);
    
    % for probable change
    % ids of uddated probable change
    ids_prob = change_prob > 0 & change_prob < 1;
    % ids for probable change before
    ids_prob_bef = change_prob < 0;
    
    % positions of probable change
    pos_prob = pos(ids_prob);
    pos_prob_bef = pos(ids_prob_bef);
    % percent probability of probable change
    pct_prob = change_prob(ids_prob);
    pct_prob_bef = -change_prob(ids_prob_bef);
    % time of probable change
    time_prob = time(ids_prob);
    time_prob_bef = time(ids_prob_bef);
    % magnitude of probable change
    mag_prob = shap_mag(:,ids_prob);
    mag_prob_bef = shap_mag(:,ids_prob_bef);
    
    % pass if there is no confirmed change
    l_pos = length(pos_break);
    if l_pos > 0  
        for i=1:l_pos
            [I,J] = ind2sub(jiDim,pos_break(i));
            ChangeNum(J,I) = ChangeNum(J,I) + 1;
%             % cannot record change time & magnitude due to matrix size
%             if ChangeNum(J,I) <= max_n 
                vec_time = datevecmx(time_break(i));
                ChangeTime(J,I,ChangeNum(J,I)) = vec_time(1);
                ChangeMag(J,I,ChangeNum(J,I)) = norm(mag_break(:,i));
%             end
        end
    end
    
    % pass if there is no probable change
    l_pos = length(pos_prob);
    if l_pos > 0
        for i=1:l_pos
            [I,J] = ind2sub(jiDim,pos_prob(i));
            PChangeProb(J,I) = uint8(100*pct_prob(i));
            vec_time = datevecmx(time_prob(i));
            PChangeTime(J,I) = vec_time(1);
            PChangeMag(J,I) = norm(mag_prob(:,i));
        end
    end
       
    % pass if there is no probable change before
    l_pos = length(pos_prob_bef);
    if l_pos > 0
        for i=1:l_pos
            [I,J] = ind2sub(jiDim,pos_prob_bef(i));
            PChangeProbBef(J,I) = uint8(100*pct_prob_bef(i));
            vec_time = datevecmx(time_prob_bef(i));
            PChangeTimeBef(J,I) = vec_time(1);
            PChangeMagBef(J,I) = norm(mag_prob_bef(:,i));
        end
    end
end

% actual maximum change number
max_n = max(ChangeNum(:));
% resize time matrix
ChangeTime = ChangeTime(:,:,1:max_n);
% resize magnitude matrix
ChangeMag = ChangeMag(:,:,1:max_n);

% write ENVI files
% confirmed change
enviwrite([v_input.l_dir,'/',n_map,'/ChangeNum'],ChangeNum,'uint8',res,jiUL,'bsq',zc);
enviwrite([v_input.l_dir,'/',n_map,'/ChangeTime'],ChangeTime,'uint16',res,jiUL,'bsq',zc);
enviwrite([v_input.l_dir,'/',n_map,'/ChangeMag'],ChangeMag,'single',res,jiUL,'bsq',zc);

% probable change after
enviwrite([v_input.l_dir,'/',n_map,'/PChangeProbAft'],PChangeProb,'uint8',res,jiUL,'bsq',zc);
enviwrite([v_input.l_dir,'/',n_map,'/PChangeTimeAft'],PChangeTime,'uint16',res,jiUL,'bsq',zc);
enviwrite([v_input.l_dir,'/',n_map,'/PChangeMagAft'],PChangeMag,'single',res,jiUL,'bsq',zc);

% probable change before
enviwrite([v_input.l_dir,'/',n_map,'/PChangeProbBef'],PChangeProbBef,'uint8',res,jiUL,'bsq',zc);
enviwrite([v_input.l_dir,'/',n_map,'/PChangeTimeBef'],PChangeTimeBef,'uint16',res,jiUL,'bsq',zc);
enviwrite([v_input.l_dir,'/',n_map,'/PChangeMagBef'],PChangeMagBef,'single',res,jiUL,'bsq',zc);




