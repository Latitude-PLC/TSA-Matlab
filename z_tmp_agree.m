% Function for comparing agreement between CCDC and Trends
%% agreement within Ecoregion
clc
pwd
ref_eco = enviread('Eco_Mosaic');
ref_ard = enviread('Eco_Mosaic_ARD');
% ref_ard = enviread('Eco_Mosaic_ARD_Refine');

idv_eco = ref_eco > 0 & ref_eco ~= 3 & ref_eco ~= 10;% & ref_ard > 0;
idv_ard = ref_ard > 0 & ref_ard ~= 3 & ref_ard ~= 10;

sum_eco = sum(idv_eco(:));
sum_ard = sum(idv_ard(:));

im_eco = enviread('Cover_Mosaic');
im_eco = im_eco(:,:,16); % 2000 map

im_ard = enviread('Cover_Mosaic_ARD_v2');
im_ardv2 = im_ard(:,:,16);

im_ard = enviread('Cover_Mosaic_ARD_v3');
im_ardv3 = im_ard(:,:,16);

im_ard = enviread('Cover_Mosaic_ARD_v4');
im_ardv4 = im_ard(:,:,16);

agree_eeco = idv_eco & ref_eco == im_eco;
agree_eardv2 = idv_eco & ref_eco == im_ardv2;
agree_eardv3 = idv_eco & ref_eco == im_ardv3;
agree_eardv4 = idv_eco & ref_eco == im_ardv4;

% agreement within ecoregion
sum(agree_eeco(:))/sum_eco
sum(agree_eardv2(:))/sum_eco
sum(agree_eardv3(:))/sum_eco
sum(agree_eardv4(:))/sum_eco

%%
pwd
agree_eco = idv_ard & ref_ard == im_eco;
agree_ardv2 = idv_ard & ref_ard == im_ardv2;
agree_ardv3 = idv_ard & ref_ard == im_ardv3;
agree_ardv4 = idv_ard & ref_ard == im_ardv4;

% agreement within ecoregion
sum(agree_eco(:))/sum_ard
sum(agree_ardv2(:))/sum_ard
sum(agree_ardv3(:))/sum_ard
sum(agree_ardv4(:))/sum_ard

