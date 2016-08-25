% renaming all the files
imf=dir('L*'); % folder names
num_l=size(imf,1); % num of folder to be processed
% add '_MTLstack' to all files
add_n = '_MTLstack';
for i=1:num_l
    fprintf('Renaming the %dth folder\n',i);
    cd(imf(i).name);
    % rename all file
    
    % find stack images
    n_Fmask = dir('L*.hdr');
    % name of image
    n_Fmask_hdr = n_Fmask.name;
    % length of name
    l_n = length(n_Fmask_hdr);
    
    if l_n == 25 % No stack name 
        % add stack name to the image folder
        n_in_img = n_Fmask_hdr(1:end-4);
        n_in_hdr = n_Fmask_hdr;
        n_out_img = [n_in_img,add_n];
        n_out_hdr = [n_in_img,add_n,'.hdr'];
    elseif l_n ~= 34 % has the wrong stack name
        % revise names
        n_Fmask = dir('L*stack');
        n_Fmask_img = n_Fmask.name;
        
        n_in_img = n_Fmask_img;
        n_in_hdr = n_Fmask_hdr;
        n_out_img = [n_Fmask_img(1:21),add_n];
        n_out_hdr = [n_Fmask_img(1:21),add_n,'.hdr'];
    else
        fprintf('Length NOT right!!!\n');
    end
    
    if length(n_in_hdr) ~= length(n_out_hdr)
        movefile(n_in_hdr,n_out_hdr);
    end
    
    if length(n_in_img) ~= length(n_out_img)
        movefile(n_in_img,n_out_img);
    end
    
    fprintf('%s image \n',n_Fmask_hdr);
    cd ..
end