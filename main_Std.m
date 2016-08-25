pwd
imf=dir('L*'); % folder names
num_l=size(imf,1); % num of folder to be processed
nbs = [3,4,5];

parfor i = 1:num_l
    fprintf('processing the %d/%dth folder ...\n',i,num_l);
    
    cd(imf(i).name)
    file = dir('L*stack');
    [im,jiDim,jiUL,resolu,ZC] = enviread(file.name);

    imstd = zeros(jiDim(2),jiDim(1),length(nbs));
    for j = 1:length(nbs)
        % imstd(:,:,j) = stdfilt(double(im(:,:,nbs(j))),ones(9));
        imstd(:,:,j) = conv2(double(im(:,:,nbs(j))),ones(7)/7^2,'same');
    end
    
    % write std image to each image folder
    enviwrite('std345',imstd,'uint16',resolu,jiUL,'bip',ZC);
    
    cd ..
end