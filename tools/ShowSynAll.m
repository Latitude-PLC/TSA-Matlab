function ShowSynAll(yy,mm,dd)
% This function is used to provde synthetic data for all pixels
% two digits for QA band
% Scripts for runing jobs

% 1. inputs: the time for predicted synthetic data
% yy = 2000;
% mm = 12;
% dd = 1;
% OR
% yy = 2000;
% doy = 120;

% addpath('/usr3/graduate/zhuzhe/Algorithms/CCDC/Script');
v_input=main_Inputs;

if nargin == 2
    doy = mm;
    j_date = datenum(yy,1,0) + doy;
    s_date = yy*1000 + doy;
elseif nargin == 3
    % converted to julian date
    j_date=datenum(yy,mm,dd);
    % j0_date=datenum(yy,1,0);
    % ddd=j_date-j0_date;
    % date for show i.e. 19990815
    s_date=yy*10000+mm*100+dd;
else
    fprintf('Check input variable number!\n');
    return;
end

nbands=v_input.nbands-1;% number of bands in the image
ncoefs=v_input.num_c;% number of coefficients

% dimension and projection of the image
cd(v_input.l_dir);
nrows = v_input.ijdim(1);
ncols = v_input.ijdim(2);
jiDim = [ncols,nrows];
jiUL = v_input.jiul;
res = v_input.resolu;
zc = v_input.zc;

% produce synthetic data (1, 2, 3, 4, 5, 7, and QA)
SR=-9999*ones(nrows,ncols,nbands,'int16'); % 1. change here for excluding thermal

% folder saved everything
n_rst = v_input.name_rst;
cd(n_rst);% TSFitMap

% make Predict folder for storing predict images
n_pre = v_input.name_pre;%'PredictAll';

if isempty(dir(n_pre))
    mkdir(n_pre);
end

% folder saved the recorded structures
n_rec = v_input.name_rec;
cd(n_rec);

imf=dir('record_change*'); % folder names
num_line=size(imf,1);

for line=1:num_line
    fprintf('Processing %.2f percent\n',100*(line/num_line));
    load(imf(line).name);
    
    % postions & coefficients
    pos=[rec_cg.pos];
    % continue if there is no model available
    l_pos=length(pos);
    if l_pos==0
        continue;
    end
    % get coefficients matrix
    coefs=[rec_cg.coefs];
    % reshape coefs
    coefs=reshape(coefs,nbands*ncoefs,l_pos);
    % get category matrix
    category = [rec_cg.category];
    
    % model start, end, & break time for prediction
    model_start = [rec_cg.t_start];
    model_end = [rec_cg.t_end];
    model_break = [rec_cg.t_break];
    % model on the right
    ids_right = model_start > j_date;
    % model on the left
    ids_left = (model_end < j_date & model_break == 0)|(model_break <= j_date & model_break > 0);
    % id within model interval
    ids_middle = model_start <= j_date & ((model_end >= j_date & model_break == 0) | (model_break > j_date & model_break > 0));

    % position for model in the middle
    pos_middle = pos(ids_middle);
    % coefficients for model in the middle
    coefs_middle = coefs(:,ids_middle);
    % category for model in the middle
    category_middle = category(ids_middle);
    
    % positions for the nearest model on the right
    pos_right = pos(ids_right);
    [pos_near_right,ids_near_right] = unique(pos_right,'first');
    % coefficients for the nearest model on the right
    coefs_right = coefs(:,ids_right);
    coefs_near_right = coefs_right(:,ids_near_right);
    % category for the nearest model on the right
    category_right = category(ids_right);
    category_near_right = category_right(ids_near_right);
    
    % postions for the nearest model on the left
    pos_left = pos(ids_left);
    [pos_near_left,ids_near_left] = unique(pos_left,'last');
    % coefficients for the nearest model on the left
    coefs_left = coefs(:,ids_left);
    coefs_near_left = coefs_left(:,ids_near_left);
    % category for the nearest model on the left
    category_left = category(ids_left);
    category_near_left = category_left(ids_near_left);
    
    % pass if there is no nearest model on the left 
    l_pos=length(pos_near_left);
    if l_pos > 0  
        % providing predictions
        for i=1:l_pos
            [I,J]=ind2sub(jiDim,pos_near_left(i));
            for j_b=1:nbands-1 % excluding thermal band
                SR(J,I,j_b)=autoTSPred(j_date,coefs_near_left(((j_b-1)*ncoefs+1):j_b*ncoefs,i));
            end
            % QA band
            SR(J,I,end) = 20 + category_near_left(i); % model forward predicted values         
        end
    end
    
    % pass if there is no nearest model on the right 
    l_pos=length(pos_near_right);
    if l_pos > 0  
        % providing predictions
        for i=1:l_pos
            [I,J]=ind2sub(jiDim,pos_near_right(i));
            for j_b=1:nbands-1 % excluding thermal band
                SR(J,I,j_b)=autoTSPred(j_date,coefs_near_right(((j_b-1)*ncoefs+1):j_b*ncoefs,i));
            end
            % QA band
            SR(J,I,end) = 10 + category_near_right(i); % model backward predicted values         
        end
    end
    
    % pass if there is no nearest model in the middle 
    l_pos=length(pos_middle);
    if l_pos > 0  
        % providing estimations
        for i=1:l_pos
            [I,J]=ind2sub(jiDim,pos_middle(i));
            for j_b=1:nbands-1 % excluding thermal band
                SR(J,I,j_b)=autoTSPred(j_date,coefs_middle(((j_b-1)*ncoefs+1):j_b*ncoefs,i));
            end
            % QA band
            SR(J,I,end) = 0 + category_middle(i); % model estimated values         
        end
    end    
end

% write ENVI files
enviwrite([v_input.l_dir,'/',n_rst,'/',n_pre,'/SR',num2str(s_date)],SR,'int16',res,jiUL,'bsq',zc);






