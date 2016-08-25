%% Single-date cloud, cloud shadow, & snow masks
clear
clc

imf=dir('L*'); % folder names
num_l=size(imf,1); % num of folder to be processed
prct = zeros(num_l,1);

parfor i=1:num_l
    fprintf('Processing the %dth folder\n',i);
    cd(imf(i).name);
    n_Fmask = dir('L*stack');
    tmp_Fmask = enviread(n_Fmask.name);
    im_Fmask = tmp_Fmask(:,:,end);
    
    clear_stat = im_Fmask < 2;
    
    prct(i) = 100*sum(clear_stat(:))/length(clear_stat(:));
    
    fprintf('%.2f percent of clear observations for %s image \n',prct(i),n_Fmask.name);
    cd ..
end

save('prct','prct');