function ShowOneMap(all_yrs,mm,dd)
% This function is used to provde all change maps for each year
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

% produce confirmed change map
OneMap = zeros(nrows,ncols,max_n,'uint16'); % DOY (0~356)

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
 
    % one variable
    texture = [rec_cg.texture];
    
    % For Condition & Cover Map
    for i = 1:l_pos
        % get row and col
        [I,J] = ind2sub(jiDim,pos(i));
                                    
        % year (band) the curve belongs to
        n_band = jul_d >= t_start(i) & (jul_d <= t_end(i) | jul_d < t_break(i));
        
        % write number of clear observations to NumberMap
        OneMap(J,I,n_band) = texture(3*(i-1)+1);               
    end
end
cd ..

enviwrite([v_input.l_dir,'/',n_map,'/OneMap'],OneMap,'uint16',res,jiUL,'bip',zc);