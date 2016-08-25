% function rmse_stat
% This function is used to provde synthetic data for all pixels
% Scripts for runing jobs
addpath('~/ccdc');
v_input=main_Inputs;
% Inputs:
% image with DN values show what land cover class
clear;
v_input = ccdc_Inputs;

l_dir = v_input.l_dir;
% Constants:
% number of coefficients
num_c = v_input.num_c;
% number of bands
nbands = v_input.nbands; 
% image dimension
nrows = v_input.ijdim(1);
ncols = v_input.ijdim(2);
jiDim = [ncols,nrows];


% folder saved everything
n_rst = v_input.name_rst;
cd(n_rst);% TSFitMap

imf=dir('record_change*'); % folder names
num_line=size(imf,1);

% define rmse and mean
f_rmse = zeros(nbands-1,10^8);
f_mean = f_rmse;
j_tot = 0;
for line = 1:num_line
    fprintf('Processing %.2f percent\n',100*(line/num_line));
    load(imf(line).name);
    
    % postions & coefficients
    pos=[rec_cg.pos];
    % continue if there is no model available
    l_pos=length(pos);
    if l_pos==0
        continue;
    end

    % get RMSE
    rmse=[rec_cg.rmse];
    % coefficients
    coefs = [rec_cg.coefs];
    % start time
    t_start = [rec_cg.t_start];
    % end time
    t_end = [rec_cg.t_end];
    % reshape coefs
    coefs = reshape(coefs,num_c,nbands-1,[]);
    
    for j = 1:l_pos
        f_mean(:,j+j_tot) = coefs(1,:,j) + (t_start(j)+t_end(j))*coefs(2,:,j)/2;
        f_rmse(:,j+j_tot) = rmse(:,j);
    end
    j_tot = j_tot + l_pos;
end
f_mean = f_mean(:,1:j_tot);
f_rmse = f_rmse(:,1:j_tot);
%%
ids = randperm(size(f_mean,2));

for i_b = 1:nbands-1
   figure;
   scatter(f_mean(i_b,ids(1:1000)),f_rmse(i_b,ids(1:1000)));
end

cd ..
%%
save('f_mean','f_mean','-v7.3');
save('f_rmse','f_rmse','-v7.3');

%%
fix_itv = 100;
all_itv = (0:fix_itv:2000);
l_itv = length(all_itv);
stat_mean = zeros(nbands-1,l_itv);
stat_sum = stat_mean;

for i_b = 2:nbands-2
    i_b
    for i_itv = 1:l_itv
        ids = f_mean(i_b,:) > all_itv(i_itv) & f_mean(i_b,:) <= all_itv(i_itv) + fix_itv;
        stat_mean(i_b,i_itv) = mean(f_rmse(i_b,ids));
        stat_sum(i_b,i_itv) = sum(ids);
    end
end

%%
stat_mean = zeros(nbands-1,1);
stat_sum = stat_mean;

for i_b = 2:nbands-2
    i_b
    ids = f_mean(i_b,:) > 0 & f_mean(i_b,:) <= 500;
    stat_mean(i_b) = mean(f_rmse(i_b,ids));
    stat_sum(i_b) = sum(ids);
end

%%
for i_b = 2:nbands-2
    figure;
    plot(all_itv+fix_itv/2,stat_mean(i_b,:)./(all_itv+fix_itv/2),'r-');
    hold on
    plot(all_itv+fix_itv/2,stat_mean(i_b,:)/1000,'g-');
    legend('Relative RMSE','Absolute RMSE');
    hold on;
    hb=bar(all_itv+fix_itv/2,stat_sum(i_b,:)/max(stat_sum(i_b,:)));
    ch=get(hb,'child');
    set(ch,'facea',0.5);
end
pause;
close all;
        





