% function ShowCMap(all_yrs,mm,dd)
% This function is used to provde change change maps for each year
% Results for LCMAP

% Version 1.00 No disturbance class in the cover map (11/06/2015)
% Tools
addpath('~/ccdc');
% INPUTS:
all_yrs = 1985:2014; % all of years for producing maps
mm = 7;
dd = 1;
%
v_input = ccdc_Inputs;
pwd

% dimension and projection of the image
nrows = v_input.ijdim(1);
ncols = v_input.ijdim(2);
jiDim = [ncols,nrows];
jiUL = v_input.jiul;
res = v_input.resolu;
zc = v_input.zc;
% number of coefficients
num_c = v_input.num_c;
% number of bands
nbands = v_input.nbands;
% max number of maps
max_n = length(all_yrs);
% all julidan dates
jul_d = datenummx(all_yrs,mm,dd);

% produce confirmed change map
ChangeMap = 9999*ones(nrows,ncols,max_n,'uint16'); % DOY (0~356)
% produce confirmed change vector magnitude map
ChangeMagMap = 9999*ones(nrows,ncols,max_n,'single'); % DOY (0~inf)

% % produce land cover map
% CoverMap = 255*ones(nrows,ncols,max_n,'uint8'); % Trends categories (0~11)
% % produce land cover map QA (unsupervised emsemble margin)
% CoverQAMap = 255*ones(nrows,ncols,max_n,'uint8'); % Trends categories (0~100)
% 
% produce condition map at all_yrs mm dd
ConditionMap = 9999*ones(nrows,ncols,max_n,'single'); % EVI slope
% produce QA band for model estimation
QAMap = 255*ones(nrows,ncols,max_n,'uint8'); % QA map

% make Predict folder for storing predict images
n_map=v_input.name_map;% 'CCDCMap';
if isempty(dir(n_map))
    mkdir(n_map);
end

% cd to the folder for storing recored structure
cd(v_input.name_rst);

imf = dir('record_change*'); % folder names
num_line = size(imf,1);
line_pt = 0;

for line = 1:num_line
    if 100*(line/num_line) - line_pt > 1
        fprintf('Processing %.0f percent\n',ceil(100*(line/num_line)));
        line_pt = 100*(line/num_line);
    end
    
    % load one line of time series models
    load(imf(line).name);
    
    % postions
    pos = [rec_cg.pos];
  
    % continue if there is no model available
    l_pos = length(pos);
    if l_pos == 0
        continue
    end
    
    % matrix of each component
    % start time
    t_start = [rec_cg.t_start];
    % end time
    t_end = [rec_cg.t_end];
    % break time
    t_break = [rec_cg.t_break];
    % coefs
    coefs = [rec_cg.coefs];
    % reshape coefs
    coefs = reshape(coefs,num_c,nbands-1,[]);
    % change probability
    change_prob = [rec_cg.change_prob];
%     % class
%     class = [rec_cg.class];
    % regression QA
    categ = [rec_cg.category];
%     % classification QA
%     class_qa = [rec_cg.classQA];
    % change vector
    mag = [rec_cg.magnitude];
    mag = reshape(mag,nbands-1,[]);
    
    % For Condition & Cover Map
    for i = 1:l_pos
        % get row and col
        [I,J] = ind2sub(jiDim,pos(i));
        
        % initialize pixels have at least one model
        if sum(ChangeMap(J,I,:) == 9999) == max_n
            % write EVI slope to ConditionMap
            ConditionMap(J,I,:) = 0;
            % write QAMap
            QAMap(J,I,:) = 0;
%             % write land cover to CoverMap
%             CoverMap(J,I,:) = 0;
%             % write land cover to CoverQAMap
%             CoverQAMap(J,I,:) = 0;
            % write doy to ChangeMap
            ChangeMap(J,I,:) = 0;
            % write doy to ChangeMagMap
            ChangeMagMap(J,I,:) = 0;
        end
        
%         % give next class for disturbed class
%         if i > 1
%             if pos(i) == pos(i-1) % same location
%                 n_dist = jul_d < t_start(i) & jul_d >= t_break(i-1);
%                 % write land cover (distubed) to CoverMap
%                 CoverMap(J,I,n_dist) = class(i);
%             end
%         end
                            
        % year (band) the curve belongs to
        n_band = jul_d >= t_start(i) & (jul_d <= t_end(i) | jul_d < t_break(i));
        
        % year (band) the distubed period
        
        % get EVI values
        % start values
        b_start = coefs(1,1,i) + t_start(i)*coefs(2,1,i);
        r_start = coefs(1,3,i) + t_start(i)*coefs(2,3,i);
        n_start = coefs(1,4,i) + t_start(i)*coefs(2,4,i);
        % end values
        b_end = coefs(1,1,i) + t_end(i)*coefs(2,1,i);
        r_end = coefs(1,3,i) + t_end(i)*coefs(2,3,i);
        n_end = coefs(1,4,i) + t_end(i)*coefs(2,4,i);
        % start of EVI
        EVI_start = 2.5*(n_start - r_start)/(n_start + 6*r_start - 7.5*b_start + 10000);
        % end of EVI
        EVI_end = 2.5*(n_end - r_end)/(n_end + 6*r_end - 7.5*b_end + 10000);
        
        % EVI slope
        EVI_slope = 10000*(EVI_end - EVI_start)/(t_end(i)-t_start(i));
        
        % write EVI slope to ConditionMap
        ConditionMap(J,I,n_band) = EVI_slope;
        % write QAMap
        QAMap(J,I,n_band) = categ(i);
%         
%         % write land cover to CoverMap
%         CoverMap(J,I,n_band) = class(i);       
%         % write land cover to CoverQAMap
%         CoverQAMap(J,I,n_band) = class_qa(i);
        
        if change_prob(i) == 1
            % get the time
            vec_time = datevecmx(t_break(i));
            % get the year
            n_year = vec_time(1);
            % get day-of-year (doy)
            n_doy = t_break(i) - datenummx(n_year,1,0);
            % get the band number
            n_band = all_yrs == n_year;
            
            % write doy to ChangeMap
            ChangeMap(J,I,n_band) = n_doy;
            % write doy to ChangeMagMap
            ChangeMagMap(J,I,n_band) = norm(mag(:,i));
        end        
    end
end
cd ..

% Change Map
ARD_enviwrite_bands([v_input.l_dir,'/',n_map,'/ChangeMap'],ChangeMap,'uint16','bsq',all_yrs,'example_img');
clear ChangeMap
ARD_enviwrite_bands([v_input.l_dir,'/',n_map,'/ChangeMagMap'],ChangeMagMap,'single','bsq',all_yrs,'example_img');
clear ChangeMagMap

% Condition Map
ARD_enviwrite_bands([v_input.l_dir,'/',n_map,'/ConditionMap'],ConditionMap,'single','bsq',all_yrs,'example_img');
clear ConditionMap
ARD_enviwrite_bands([v_input.l_dir,'/',n_map,'/QAMap'],QAMap,'uint8','bsq',all_yrs,'example_img');
clear QAMap


% % write ENVI files
% 
% % Change Map
% enviwrite_bands([v_input.l_dir,'/',n_map,'/ChangeMap'],ChangeMap,'uint16',res,jiUL,'bsq',zc,all_yrs);
% clear ChangeMap
% enviwrite_bands([v_input.l_dir,'/',n_map,'/ChangeMagMap'],ChangeMagMap,'single',res,jiUL,'bsq',zc,all_yrs);
% clear ChangeMagMap
% 
% % Condition Map
% enviwrite_bands([v_input.l_dir,'/',n_map,'/ConditionMap'],ConditionMap,'single',res,jiUL,'bsq',zc,all_yrs);
% clear ConditionMap
% enviwrite_bands([v_input.l_dir,'/',n_map,'/QAMap'],QAMap,'uint8',res,jiUL,'bsq',zc,all_yrs);
% clear QAMap