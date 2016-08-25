function Fmask_stats

% Revisions: 
% v1.0 Only output the three useful Fmask statistics (03/21/2016)


v_input = ccdc_Inputs;

% dimension and projection of the image
nrows = v_input.ijdim(1);
ncols = v_input.ijdim(2);
jiDim = [ncols,nrows];
jiUL = v_input.jiul;
res = v_input.resolu;
zc = v_input.zc;

imf=dir('L*'); % folder names
num_l=size(imf,1); % num of folder to be processed

% statistics 
Fmask_stat = zeros(nrows,ncols,4);
All_stat = zeros(nrows,ncols);

% % add '_MTLstack' to all files
% add_n = '_MTLstack';
% for i=1:num_l
%     fprintf('Renaming the %dth folder\n',i);
%     cd(imf(i).name);
%     % rename all file
%     
%     % find stack images
%     n_Fmask = dir('L*.hdr');
%     % name of image
%     n_Fmask_hdr = n_Fmask.name;
%     % length of name
%     l_n = length(n_Fmask_hdr);
%     
%     if l_n == 25 % No stack name 
%         % add stack name to the image folder
%         n_in_img = n_Fmask_hdr(1:end-4);
%         n_in_hdr = n_Fmask_hdr;
%         n_out_img = [n_in_img,add_n];
%         n_out_hdr = [n_in_img,add_n,'.hdr'];
%     elseif l_n == 34 || l_n == 43 % has >= 1 stack name
%         % revise names
%         n_Fmask = dir('L*_MTLstack');
%         n_Fmask_img = n_Fmask.name;
%         
%         n_in_img = n_Fmask_img;
%         n_in_hdr = n_Fmask_hdr;
%         n_out_img = n_Fmask_img(1:30);
%         n_out_hdr = [n_Fmask_img(1:30),'.hdr'];
%     else
%         fprintf('Length NOT right!!!\n');
%     end
%     
%     if length(n_in_hdr) ~= length(n_out_hdr)
%         movefile(n_in_hdr,n_out_hdr);
%     end
%     
%     if length(n_in_img) ~= length(n_out_img)
%         movefile(n_in_img,n_out_img);
%     end
%     
%     fprintf('%s image \n',n_Fmask_hdr);
%     cd ..
% end

% go through all images
for i=1:num_l
    fprintf('Processing the %dth folder\n',i);
    cd(imf(i).name);
    n_Fmask = dir('L*stack');
    tmp_Fmask = enviread(n_Fmask.name);
    im_Fmask = tmp_Fmask(:,:,end);
    
    Fmask_stat(:,:,4) = Fmask_stat(:,:,4) + double(im_Fmask == 4);
    Fmask_stat(:,:,3) = Fmask_stat(:,:,3) + double(im_Fmask == 3);
    Fmask_stat(:,:,2) = Fmask_stat(:,:,2) + double(im_Fmask == 1);
    Fmask_stat(:,:,1) = Fmask_stat(:,:,1) + double(im_Fmask == 0);
    All_stat = All_stat + double(im_Fmask < 255);
    
    fprintf('%s image \n',n_Fmask.name);
    cd ..
end

% write into the same file
% cloud probability
Fmask_stat(:,:,4) = 100*Fmask_stat(:,:,4)./All_stat;
% snow probability
Fmask_stat(:,:,3) = 100*Fmask_stat(:,:,3)./(Fmask_stat(:,:,1)+Fmask_stat(:,:,2)+Fmask_stat(:,:,3)+0.01);
% clear water probability
Fmask_stat(:,:,2) = 100*Fmask_stat(:,:,2)./(Fmask_stat(:,:,1)+Fmask_stat(:,:,2)+0.01);
% clear land probability
Fmask_stat(:,:,1) = 100*Fmask_stat(:,:,1)./(Fmask_stat(:,:,1)+Fmask_stat(:,:,2)+0.01);

Fmask_stat = uint8(Fmask_stat(:,:,[2,3,4]));

% write ENVI files
% enviwrite([v_input.l_dir,'/ANC/Fmask_stat'],Fmask_stat,'uint8',res,jiUL,'bip',zc);
ARD_enviwrite([v_input.l_dir,'/ANC/Fmask_stat'],Fmask_stat,'uint8','bip','example_img');
end