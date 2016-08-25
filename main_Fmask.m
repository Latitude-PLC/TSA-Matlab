pwd
outimf=dir('newARD*'); % folder names
outnum_l=size(outimf,1); % num of folder to be processed

parfor i = 2:outnum_l
    cd(outimf(i).name)
    Fmask_stats;
    cd ..
end