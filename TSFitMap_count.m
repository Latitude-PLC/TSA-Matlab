% function main_Fmask(location_n)

imf=dir('newARD*'); % folder names
num_l=size(imf,1); % num of folder to be processed

for i=1:num_l
    cd(imf(i).name);
    cd TSFitMap
    
    n = dir('record*');
    n = size(n,1);
    
    if n < 5000
        n
        pwd
    end

    cd ..
    cd ..;
end