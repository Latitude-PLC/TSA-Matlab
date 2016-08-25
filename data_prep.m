function data_prep(task,ntask)
% This code prepares the espa Landsat data into CCDC format
% Able to handle both TIF and ENVI format
% Selecting the surface reflectance and temperataure, Fmask results and
% stacking them into bip and save the results into the images folder
%
% Version 1.0 data preparation code for Landsat data from espa (Zhe Zhu 02/05/2015)
%
% Input: 
% task: the number of core to be used
% ntask: total number of cores to be used
%
% Output:
% data prepared for CCDC
%
% Example:
% 1. cd to the folder where all tar.zip files are stored
% e.g. cd /projectnb/landsat/projects/LCMAP/p023r037
% 2. if using three cores in each Matlab command window
% data_prep(1,3);
% data_prep(2,3);
% data_prep(3,3);
%
% print to screen
pwd
% task = 1;
% ntask = 1;
%
% current directory
dir_cur = pwd; 
% get all zip file names
dir_f = dir('L*.tar.gz');
% number of files
n_f = size(dir_f,1);
% total number of bands
n_bs = 8;

% divid jobs
ids_all = mod(1:n_f,ntask);
ids_all(ids_all == 0) = ntask;
% find ids that are divided by ntask and equal to task
ids_div = find(ids_all == task);
% name of the temporary folder for extracting zip files
name_tmp = 'tmp';
n_tmp = ['/',name_tmp,num2str(task),'/'];

for i = ids_div
    
    % check if folder exsit or not
    % names of image folder that are processed
    n_img = dir([dir_cur,'/images/','L*']);
    num_img = size(n_img,1);
    % check all folders we have
    % record exist or not
    rec_exist = 0;
    for i_check = 1:num_img
        % each image folder name
        tmp_img = n_img(i_check).name;        
        % the ith zip file
        tmp_zip = dir_f(i).name;        
        if strcmp(tmp_img(1:16),tmp_zip(1:16)) == true
            rec_exist = 1;
            break;
        end
    end
    % continue if the folder already exist
    if rec_exist > 0
        fprintf('%s exixt in images folder\n',tmp_img);
        continue;
    end
    
    if str2num(dir_f(i).name(3)) < 8 % Landsats 4-7
        
        % unzip the images to a temporary folder
        fprintf('Unzip the %dth image ...\n',i);
        n_gun = gunzip(dir_f(i).name,[name_tmp,num2str(task)]);
        n_tar = untar([dir_cur,'/',char(n_gun)],[name_tmp,num2str(task)]);
        
        % decide image format (tif or envi)
        env_cfmask = dir([dir_cur,n_tmp,'L*cfmask.img']);
        tif_cfmask = dir([dir_cur,n_tmp,'L*cfmask.tif']);
        
        if ~isempty(env_cfmask) % envi format
            % picking surf ref 1-7, bt, and cfmask and save to the images folder
            % get names of surf 1-7
            fprintf('Reading images ...\n');
            
            % read cfmask first to caculate clear pixel percet
            env_cfmask = [dir_cur,n_tmp,env_cfmask.name];
            [cfmask,jidim,jiul,resolu,zc] = enviread(env_cfmask);
            clr_pct = sum(cfmask(:)<=1)/sum(cfmask(:)<255);
            
            if clr_pct < 0.2 % less than 20% clear observations
                % remove the tmp folder
                fprintf('Clear observation less than 20 percent (%.2f) ...\n',clr_pct*100);
                rmdir([name_tmp,num2str(task)],'s');
                continue;
            else
                % prelocate image for the stacked image
                stack = zeros(jidim(2),jidim(1),n_bs,'int16');
                stack(:,:,end) = cfmask;
            end
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band1.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b1 = enviread(n_surf);
            stack(:,:,1) = surf_b1;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band2.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b2 = enviread(n_surf);
            stack(:,:,2) = surf_b2;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band3.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b3 = enviread(n_surf);
            stack(:,:,3) = surf_b3;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band4.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b4 = enviread(n_surf);
            stack(:,:,4) = surf_b4;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band5.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b5 = enviread(n_surf);
            stack(:,:,5) = surf_b5;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band7.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b7 = enviread(n_surf);
            stack(:,:,6) = surf_b7;
            
            n_surf = dir([dir_cur,n_tmp,'L*toa_band6.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b6 = enviread(n_surf);
            stack(:,:,7) = surf_b6;
            
        elseif ~isempty(tif_cfmask) % tif format
            % picking surf ref 1-7, bt, and cfmask and save to the images folder
            % get names of surf 1-7
            fprintf('Reading images ...\n');
            
            % read cfmask first to caculate clear pixel percet
            tif_cfmask = [dir_cur,n_tmp,tif_cfmask.name];
            cfmask = geotiffread(tif_cfmask);
            clr_pct = sum(cfmask(:)<=1)/sum(cfmask(:)<255);
            
            if clr_pct < 0.2 % less than 20% clear observations
                % remove the tmp folder
                fprintf('Clear observation less than 20 percent (%.2f) ...\n',clr_pct*100);
                rmdir([name_tmp,num2str(task)],'s');
                continue;
            else
                % get projection information from geotiffinfo
                info = geotiffinfo(tif_cfmask);
                jidim = [info.SpatialRef.RasterSize(2),info.SpatialRef.RasterSize(1)];
                jiul = [info.SpatialRef.XLimWorld(1),info.SpatialRef.YLimWorld(2)];
                resolu = [info.PixelScale(1),info.PixelScale(2)];
                zc = info.Zone;
                
                % prelocate image for the stacked image
                stack = zeros(jidim(2),jidim(1),n_bs,'int16');
                % give cfmask to the last band
                stack(:,:,end) = cfmask;
            end
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band1.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b1 = geotiffread(n_surf);
            stack(:,:,1) = surf_b1;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band2.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b2 = geotiffread(n_surf);
            stack(:,:,2) = surf_b2;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band3.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b3 = geotiffread(n_surf);
            stack(:,:,3) = surf_b3;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band4.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b4 = geotiffread(n_surf);
            stack(:,:,4) = surf_b4;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band5.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b5 = geotiffread(n_surf);
            stack(:,:,5) = surf_b5;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band7.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b7 = geotiffread(n_surf);
            stack(:,:,6) = surf_b7;
            
            n_surf = dir([dir_cur,n_tmp,'L*toa_band6.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b6 = geotiffread(n_surf);
            stack(:,:,7) = surf_b6;
            
        else
            warning('There is no TIF or ENVI format cfmask!\n');
            continue;
        end
        
%         % new stacked bip image
%         % get image name from MTL file
%         n_mtl = dir([dir_cur,n_tmp,'L*MTL.txt']);
%         n_mtl = strsplit(char(n_mtl.name),'.'); % remove .txt
%         n_mtl = n_mtl(1);
%         n_stack = [char(n_mtl),'stack'];
        
        % new stacked bip image
        % get image name from cfmask file
        n_mtl = dir([dir_cur,n_tmp,'L*_cfmask.img']);
        n_mtl = strsplit(char(n_mtl.name),'_'); % remove .txt
        n_mtl = n_mtl(1);
        n_stack = [char(n_mtl),'_MTLstack'];
        
        % generate folder name
        n_mtl = strsplit(char(n_mtl),'_'); % remove _MTL
        n_mtl = char(n_mtl(1));
        
        % add directory
        n_dir = [dir_cur,'/images/',n_mtl];
        mkdir(n_dir);
        
        % write to images folder
        fprintf('Writing %s image ...\n',n_mtl);
        n_stack = [n_dir,'/',n_stack];
        enviwrite(n_stack,stack,'int16',resolu,jiul,'bip',zc);
        
        % remove the tmp folder
        rmdir([name_tmp,num2str(task)],'s');
        
    else % Landsat 8 images
        % unzip the images to a temporary folder
        fprintf('Unzip the %dth image ...\n',i);
        n_gun = gunzip(dir_f(i).name,[name_tmp,num2str(task)]);
        n_tar = untar([dir_cur,'/',char(n_gun)],[name_tmp,num2str(task)]);
        
        % decide image format tif or envi
        env_cfmask = dir([dir_cur,n_tmp,'L*cfmask.img']);
        tif_cfmask = dir([dir_cur,n_tmp,'L*cfmask.tif']);
        
        if ~isempty(env_cfmask) % envi format
            % picking surf ref 1-7, bt, and cfmask and save to the images folder
            % get names of surf 1-7
            fprintf('Reading images ...\n');
            
            % read cfmask first to caculate clear pixel percet
            env_cfmask = [dir_cur,n_tmp,env_cfmask.name];
            [cfmask,jidim,jiul,resolu,zc] = enviread(env_cfmask);
            clr_pct = sum(cfmask(:)<=1)/sum(cfmask(:)<255);
            
            if clr_pct < 0.2 % less than 20% clear observations
                % remove the tmp folder
                fprintf('Clear observation less than 20 percent (%.2f) ...\n',clr_pct*100);
                rmdir([name_tmp,num2str(task)],'s');
                continue;
            else
                % prelocate image for the stacked image
                stack = zeros(jidim(2),jidim(1),n_bs,'int16');
                stack(:,:,end) = cfmask;
            end
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band2.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b1 = enviread(n_surf);
            stack(:,:,1) = surf_b1;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band3.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b2 = enviread(n_surf);
            stack(:,:,2) = surf_b2;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band4.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b3 = enviread(n_surf);
            stack(:,:,3) = surf_b3;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band5.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b4 = enviread(n_surf);
            stack(:,:,4) = surf_b4;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band6.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b5 = enviread(n_surf);
            stack(:,:,5) = surf_b5;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band7.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b7 = enviread(n_surf);
            stack(:,:,6) = surf_b7;
            
            n_surf = dir([dir_cur,n_tmp,'L*toa_band10.img']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b6 = enviread(n_surf);
            stack(:,:,7) = surf_b6;
            
        elseif ~isempty(tif_cfmask) % tif format
            % picking surf ref 1-7, bt, and cfmask and save to the images folder
            % get names of surf 1-7
            fprintf('Reading images ...\n');
            
            % read cfmask first to caculate clear pixel percet
            tif_cfmask = [dir_cur,n_tmp,tif_cfmask.name];
            cfmask = geotiffread(tif_cfmask);
            clr_pct = sum(cfmask(:)<=1)/sum(cfmask(:)<255);
            
            if clr_pct < 0.2 % less than 20% clear observations
                % remove the tmp folder
                fprintf('Clear observation less than 20 percent (%.2f) ...\n',clr_pct*100);
                rmdir([name_tmp,num2str(task)],'s');
                continue;
            else                
                % get projection information from geotiffinfo
                info = geotiffinfo(tif_cfmask);
                jidim = [info.SpatialRef.RasterSize(2),info.SpatialRef.RasterSize(1)];
                jiul = [info.SpatialRef.XLimWorld(1),info.SpatialRef.YLimWorld(2)];
                resolu = [info.PixelScale(1),info.PixelScale(2)];
                zc = info.Zone;
                
                % prelocate image for the stacked image
                stack = zeros(jidim(2),jidim(1),n_bs,'int16');
                % give cfmask to the last band
                stack(:,:,end) = cfmask;
            end
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band2.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b1 = geotiffread(n_surf);
            stack(:,:,1) = surf_b1;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band3.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b2 = geotiffread(n_surf);
            stack(:,:,2) = surf_b2;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band4.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b3 = geotiffread(n_surf);
            stack(:,:,3) = surf_b3;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band5.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b4 = geotiffread(n_surf);
            stack(:,:,4) = surf_b4;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band6.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b5 = geotiffread(n_surf);
            stack(:,:,5) = surf_b5;
            
            n_surf = dir([dir_cur,n_tmp,'L*sr_band7.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b7 = geotiffread(n_surf);
            stack(:,:,6) = surf_b7;
            
            n_surf = dir([dir_cur,n_tmp,'L*toa_band10.tif']);
            n_surf = [dir_cur,n_tmp,n_surf.name];
            surf_b6 = geotiffread(n_surf);
            stack(:,:,7) = surf_b6;
        else
            warning('There is no TIF or ENVI format cfmask!\n');
            continue;
        end
        
%         % new stacked bip image
%         % get image name from MTL file
%         n_mtl = dir([dir_cur,n_tmp,'L*MTL.txt']);
%         n_mtl = strsplit(char(n_mtl.name),'.'); % remove .txt
%         n_mtl = n_mtl(1);
%         n_stack = [char(n_mtl),'stack'];
        
        % new stacked bip image
        % get image name from cfmask file
        n_mtl = dir([dir_cur,n_tmp,'L*_cfmask.img']);
        n_mtl = strsplit(char(n_mtl.name),'_'); % remove .txt
        n_mtl = n_mtl(1);
        n_stack = [char(n_mtl),'_MTLstack'];
        
        % generate folder name
        n_mtl = strsplit(char(n_mtl),'_'); % remove _MTL
        n_mtl = char(n_mtl(1));
        
        % add directory
        n_dir = [dir_cur,'/images/',n_mtl];
        mkdir(n_dir);
        
        % write to images folder
        fprintf('Writing %s image ...\n',n_mtl);
        n_stack = [n_dir,'/',n_stack];
        enviwrite(n_stack,stack,'int16',resolu,jiul,'bip',zc);
        
        % remove the tmp folder
        rmdir([name_tmp,num2str(task)],'s');
    end
end

if task == 1
    % generate a false reference image (last cfmask band) at the end
    n_example = [dir_cur,'/images/example_img'];
    enviwrite(n_example,cfmask,'uint8',resolu,jiul,'bsq',zc);
end
% exit
end
