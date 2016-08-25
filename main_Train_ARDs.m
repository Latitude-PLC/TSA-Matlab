% Function for automatic training for ARDs
% Prepare data for Training_Strategy.m
%
% CCDC 1.4 version - Zhe Zhu, EROS, USGS
%
% Revisions: $ Date: 11/25/2015 $ Copyright: Zhe Zhu
%
% Version 1.4  Add ancillary data from NLCD and Fmask (11/25/2015)
% Version 1.3  Select sample based on best strategy (06/30/2015)
% Version 1.2  Only use undisturbed data (05/11/2015)
% Version 1.1  Use version 7.3 for storing RF model (01/10/2015)
% Version 1.0  Fixed a bug in picking the wrong pixel for training (11/08/2014)
%
% Inputs:
% image with DN values show what land cover class (Trends)
% image with DN values show what block the data are from (Trends_ids)

%% Prepare for the inputs
clear
clc

addpath('~/ccdc');
imf=dir('newARD*'); % folder names
num_l=size(imf,1); % num of folder to be processed
% how many digit for numbering grid
n_d = 2;
nfolder = imf(1).name(1:end-n_d);

% Extracting Xs & Ys
for i = 1:num_l
    fprintf('Extracting the %s folder\n',imf(i).name);
    cd(imf(i).name);
    
    % remove all previous stored Xs & Ys
    delete('Xs_*');
    delete('Ys_*');
    
    % extracting Xs & Ys and store them with unique names
    main_Select_Train_Trends1_6_part1
    cd ..
end

pwd

% Copying the 8-connected neighborhood grids
grids = [2,3,4,5;7,8,9,10;11,12,13,14;16,17,18,19];
[nr,nc] = size(grids);

% copying the selected ARD Xs & Ys
for i = 1:num_l
    [row,col] = find(grids == grids(i));
    w_r = 1;
    w_c = 1;
    wr = 1;
    wc = 1;
    
    if row - w_r < 1
        w_r = 0;
    end
    if col - w_c < 1
        w_c = 0;
    end
    if row + wr > nr
        wr = 0;
    end
    if col + wc > nc
        wc = 0;
    end
    
    % neigbourhood grids
    nb_grids = grids(row-w_r:row+wr,col-w_c:col+wc);
    
    % number of grids used
    n_times = length(nb_grids(:));
    
    % current grid in str
    if grids(i) < 10
        cngrid = ['0',num2str(grids(i))];
    else
        cngrid = num2str(grids(i));
    end
    
    % merge Xs & Ys from neighbouring folders
    Merge_Sample(cngrid,nb_grids,nfolder);
    
    cngrid %#ok<NOPTS>
    nb_grids %#ok<NOPTS>
    
    cd([nfolder,cngrid])
    % train each folder
    main_Select_Train_Trends1_6_part2(n_times)
    cd ..
end
