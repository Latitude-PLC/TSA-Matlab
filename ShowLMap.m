function ShowLMap(all_yrs,mm,dd)
% This function is used to provde all land classification maps for each year
% Results for LCMAP

% Version 1.00 No disturbance class in the cover map (11/06/2015)

% INPUTS:
% all_yrs = 1985:2014; % all of years for producing maps
addpath('~/ccdc');
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
jul_start = datenummx(all_yrs,1,1);
jul_end = datenummx(all_yrs,12,31);

% produce land cover map
CoverMap = 255*ones(nrows,ncols,max_n,'uint8'); % Trends categories (0~11)
% produce land cover map QA (unsupervised emsemble margin)
CoverQAMap = 255*ones(nrows,ncols,max_n,'uint8'); % Trends categories (0~100)

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
    % class
    class = [rec_cg.class];
    % classification QA
    class_qa = [rec_cg.classQA];
    % regression QA
    categ = [rec_cg.category];
    
    % For Condition & Cover Map
    for i = 1:l_pos
        % get row and col
        [I,J] = ind2sub(jiDim,pos(i));
        
        % initialize pixels have at least one model
        if sum(CoverMap(J,I,:) == 255) == max_n
            % write land cover to CoverMap
            CoverMap(J,I,:) = 0;
            % write land cover to CoverQAMap
            CoverQAMap(J,I,:) = 0;
        end
        
        % give disturbed class for each year
        
        % year (band) the curve belongs to
        n_band = jul_d >= t_start(i) & (jul_d <= t_end(i) | jul_d < t_break(i));        
        % write land cover to CoverMap
        CoverMap(J,I,n_band) = class(i);
        % write land cover to CoverQAMap
        CoverQAMap(J,I,n_band) = class_qa(i);        
        
        % give next land cover category for gaps
        if i > 1
            if pos(i) == pos(i-1) % same location
                n_dist = jul_d < t_start(i) & jul_d >= t_break(i-1);
                % write land cover (distubed) to CoverMap
                CoverMap(J,I,n_dist) = class(i);
                % write land cover to CoverQAMap
                CoverQAMap(J,I,n_dist) = class_qa(i);
            end
        end
        
        % year (band) the disturbed class (3)
        d_band =  t_break(i) >= jul_start & t_break(i) <= jul_end;
        % write land cover to CoverMap
        CoverMap(J,I,d_band) = 3;
        % write land cover to CoverQAMap
        CoverQAMap(J,I,d_band) = 100;
        
%         % year (band) the curve belongs to
%         n_band = jul_d >= t_start(i) & (jul_d <= t_end(i) | jul_d < t_break(i));
%         
%         % write land cover to CoverMap
%         CoverMap(J,I,n_band) = class(i);
%         % write land cover to CoverQAMap
%         CoverQAMap(J,I,n_band) = class_qa(i);

        
%         % give next class for disturbed class (v13.07)
%         if i > 1
%             if pos(i) == pos(i-1) && categ(i-1) > 30 && categ(i-1) < 40 % same location and i-1 is disturbed
%                 n_dist = jul_d < t_start(i) & jul_d >= t_start(i-1) ;
%                 % write land cover (distubed) to CoverMap
%                 CoverMap(J,I,n_dist) = class(i);
%             end
%         end
%            
%         if categ(i) < 30 || categ(i) > 40 % do not show disturbed class
%             % year (band) the curve belongs to
%             n_band = jul_d >= t_start(i) & (jul_d <= t_end(i) | jul_d < t_break(i));
%             
%             % write land cover to CoverMap
%             CoverMap(J,I,n_band) = class(i);
%             % write land cover to CoverQAMap
%             CoverQAMap(J,I,n_band) = class_qa(i);
%         end
    end
end
cd ..

% write ENVI files for ARD
% Cover Map
ARD_enviwrite_bands([v_input.l_dir,'/',n_map,'/CoverMap'],CoverMap,'uint8','bsq',all_yrs,'example_img');
clear CoverMap
ARD_enviwrite_bands([v_input.l_dir,'/',n_map,'/CoverQAMap'],CoverQAMap,'uint8','bsq',all_yrs,'example_img');
clear CoverQAMap

% % write ENVI files
% % Cover Map
% enviwrite_bands([v_input.l_dir,'/',n_map,'/CoverMap'],CoverMap,'uint8',res,jiUL,'bsq',zc,all_yrs);
% clear CoverMap
% enviwrite_bands([v_input.l_dir,'/',n_map,'/CoverQAMap'],CoverQAMap,'uint8',res,jiUL,'bsq',zc,all_yrs);
% clear CoverQAMap



