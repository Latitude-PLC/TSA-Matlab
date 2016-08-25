% Function for merging new samples
function Merge_New_Sample(cngrid,nb_grids,nfolder) 

% total samples
n_tot = length(nb_grids(:));

% merged Xs and Ys
for i = 1:n_tot
    % load(name_Xs(i).name);
    if nb_grids(i) < 10
        ngrid = ['0',num2str(nb_grid(i))];
    else
        ngrid = num2str(nb_grid(i));
    end
    
    
    sprf = [nfolder,ngrid,'/','Xs_',nfolder,ngrid];
    
    fprintf('Neighboor Grid %s\n',sprf);
    
    load([nfolder,ngrid,'/','Xs_',nfolder,ngrid,'.mat'])
    
    % load(name_Ys(i).name);
    load([nfolder,ngrid,'/','Ys_',nfolder,ngrid,'.mat'])
    
    if i == 1
        Xs = X; %#ok<*NODEF>
        Ys = Y;
        size(Y)
    else
        Xs = [Xs;X];
        Ys = [Ys;Y];
        size(Y)
    end
end
Y = Ys;
X = Xs;
size(Y)

sprf = [nfolder,cngrid,'/Xs'];
fprintf('Central Grid %s\n',sprf);

% save Xs and Ys
save([nfolder,cngrid,'/Ys'],'Y','-v7.3');
save([nfolder,cngrid,'/Xs'],'X','-v7.3');
end
    
