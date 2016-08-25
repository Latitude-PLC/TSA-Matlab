
% This function is used to provde Change Maps for each pixels
% Revisions: $ Date: 02/10/2015 $ Copyright: Zhe Zhu, EROS USGS
% Version 1.0: Show change maps (01/10/2015)
clear
v_input = ccdc_Inputs;
v_lct = pwd;

% dimension and projection of the image
nrows = v_input.ijdim(1);
ncols = v_input.ijdim(2);
jiDim = [ncols,nrows];
jiUL = v_input.jiul;
res = v_input.resolu;
zc = v_input.zc;
nbands = v_input.nbands-1;


% cd to the folder for storing recored structure
cd(v_input.name_rst);

imf=dir('record_change*'); % folder names
num_line=size(imf,1);
line_pt = 0;
for line=1:num_line
    
    if 100*(line/num_line) - line_pt > 5
        fprintf('Processing %.0f percent\n',ceil(100*(line/num_line)));
        fprintf('%s\n',v_lct);
        line_pt = 100*(line/num_line);
    end
    
    load(imf(line).name);
    
    % postions 
    pos = [rec_cg.pos];  
    
    % rmse
    rmse = [rec_cg.rmse];
    
    if length(pos) ~= size(rmse,2)
        fprintf('pos and rmse length do not match!\n');
    else
        
        if isreal(rmse) == 1
            continue
        else
            for i_pos = 1:length(pos)
                if isreal(rmse(:,i_pos)) == 0
                    [I,J] = ind2sub(jiDim,pos(i_pos));
                    fprintf('RMSE in Col = %d & Row = %d is not real!\n',I,J);
                end
            end
        end
    end
end

