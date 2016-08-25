function ShowClassMap(yy,mm,dd)
% This function is used to provde Change Type Maps for each pixels
% Two kind of maps (original and gap-filled)
 % Map1 (t_start) Map2 (t_end) Map3 (pos) 
 % Map4 (class) Map5 (change_prob) Map6 (t_break)
 
addpath('~/ccdc');
v_input = ccdc_Inputs;
pwd
% nbands=v_input.nbands-1;% number of bands in the image (except the Fmask)
% ncoefs=v_input.num_c;% number of coefficients

% dimension and projection of the image
nrows = v_input.ijdim(1);
ncols = v_input.ijdim(2);
jiDim = [ncols,nrows];
jiUL = v_input.jiul;
res = v_input.resolu;
zc = v_input.zc;
% converted to julian date
j_date=datenum(yy,mm,dd);
% date for show i.e. 19990815
s_date=yy*10000+mm*100+dd;

% produce land cover map at any date
ClassType = zeros(nrows,ncols,1,'uint8'); 
% transpose matrix
tClassType = ClassType';

% produce gap-filled land cover map at any date
fClassType = zeros(nrows,ncols,1,'uint8'); 
% transpose matrix
tfClassType = fClassType';

% left land cover (<) transpose matrix
ltClassType = tfClassType;

% make Predict folder for classification maps 
n_map = v_input.name_map;% 'CCDCMap';
if isempty(dir(n_map))
    mkdir(n_map);
end

% folder saved database
cd(v_input.name_rst);% TSFitMap
imf=dir('TSFitMap*'); % folder names
num_line=size(imf,1);

for line=1:num_line
    fprintf('Processing %.2f percent\n',100*(line/num_line));
    load(imf(line).name);
    
    % ids within the time range
    ids = Map(1,:) <= j_date & Map(2,:) >= j_date;   %#ok<NODEF>
    tClassType(Map(3,ids)) = Map(4,ids);
    
    % ids within or later than the time range
    ids = Map(2,:) >= j_date;  
    pos_right = Map(3,ids);
    class_right = Map(4,ids);
    
    [pos_near_right,ids_near_right] = unique(pos_right,'first');
    tfClassType(pos_near_right) = class_right(ids_near_right);
    
    % ids before the time range
    ids = Map(2,:) < j_date;
    pos_left = Map(3,ids);
    class_left = Map(4,ids);
    
    [pos_near_left,ids_near_left] = unique(pos_left,'last');
    ltClassType(pos_near_left) = class_left(ids_near_left);   
end

% Gap-fill priority 1) within model 2) later model 3) previous model
tfClassType(tfClassType == 0) = ltClassType(tfClassType == 0);

% tranpose back
ClassType = tClassType';
fClassType = tfClassType';

% write ENVI files
enviwrite([v_input.l_dir,'/',n_map,'/Trends_',...
    num2str(s_date)],ClassType,'uint8',res,jiUL,'bsq',zc);
enviwrite([v_input.l_dir,'/',n_map,'/Trends_GapFill_',...
    num2str(s_date)],fClassType,'uint8',res,jiUL,'bsq',zc);