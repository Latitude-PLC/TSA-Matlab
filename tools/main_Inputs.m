function v_input=main_Inputs
% function v_input=Variable_Inputs
% Function for input variables and paths
% Use the default fmask toobox developed by Zhe Zhu
addpath('~/ccdc');
% Tools of RFC
addpath('~/Algorithms/CCDC/Tools/RF_Class_C');
%% Inputs:
% 1. Location of "im_roi" & place to save "modelRF"
% Canada beetle infestation project
% v_input.l_dir='/projectnb/landsat/projects/Beetle_Infest_Map/images/';
% LCMS project
v_input.l_dir='/projectnb/landsat/projects/LCMS/4530/images/'; % done 
% v_input.l_dir='/projectnb/landsat/projects/LCMS/2727/images/'; % done 13
% v_input.l_dir='/projectnb/landsat/projects/LCMS/1228/images/'; % done
% v_input.l_dir='/projectnb/landsat/projects/LCMS/1432/images/'; % done
% v_input.l_dir='/projectnb/landsat/projects/LCMS/1637/images/'; % done
% v_input.l_dir='/projectnb/landsat/projects/LCMS/3532/images/'; % done 9

% Guangzhou project
% v_input.l_dir='/projectnb/landsat/users/fuyc/p122r044/images/'; 

% ACRE project
% v_input.l_dir='/projectnb/landsat/projects/ACRE/stacks/p004r066/images/';
% v_input.l_dir='/projectnb/landsat/projects/ACRE/stacks/p005r066/images/';
% v_input.l_dir='/projectnb/landsat/projects/ACRE/stacks/p003r067/images/';
% v_input.l_dir='/projectnb/landsat/projects/ACRE/stacks/p002r067/images/';

% 2. ground truth time interval for continuous classification
% v_input.gt=[datenum(2005,1,0),datenum(2007,12,31)]; % Boston
% v_input.gt=[datenum(2001,4,17),datenum(2001,10,26)]; % Amazon
% v_input.gt=[datenum(2000,1,1),datenum(2000,12,31)]; % LCMS
% v_input.gt=[datenum(2000,1,1),datenum(2000,12,31)]; % Composite yearly
% v_input.gt=[datenum(2000,6,1),datenum(2000,9,1)]; % Composite seasonly (June~August)
% v_input.gt=[datenum(2008,1,1),datenum(2009,12,31)]; % Guangzhou
v_input.gt=[datenum(2006,1,1),datenum(2006,12,31)]; % nlcd 3532

%% Constants:
% number of maximum coefficients
v_input.num_c = 8; 
% name of ground truth land cover map
v_input.name_roi = 'example_img'; 
% v_input.name_roi = 'training_data'; 
% folder name of all CCDC resultls 
v_input.name_rst = 'TSFitMap';
% folder name of CCDC maps
v_input.name_map = 'CCDCMap';
% folder name of CCDC predicted surf ref
v_input.name_pre = 'PredictAll';
% v_input.name_roi='class_map'; 
% get to the data folder
cd(v_input.l_dir);
% read in image dimension from reference image
[jiDim,jiul,resolu,zc] = envihdrread(v_input.name_roi);
% dimension of image [row,col]
v_input.ijdim = [jiDim(2),jiDim(1)];
% upper left coordinates
v_input.jiul = jiul;
% resolution of the image
v_input.resolu = resolu;
% zone code
v_input.zc = zc;
% number of bands (including the last Fmask band)
v_input.nbands = 8;
end