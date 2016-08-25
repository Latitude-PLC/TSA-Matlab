pwd
outimf=dir('newARD*'); % folder names
outnum_l=size(outimf,1); % num of folder to be processed

for i = 1:outnum_l
    cd(outimf(i).name)
    % ShowAccuChangeMap; % accumulated change maps
    ShowChangeMap; % annual change maps
    % ShowLMap(1985:2014,7,1); % annual land cover maps
    cd ..
end