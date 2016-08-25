function main_AdjustRMSE

% Get variable and path 
dir_l = pwd; 

% get num of total folders start with "L"
imf=dir('L*'); % folder names
% number of images
num_imgs = size(imf,1);
% filter for Landsat folders
imf = regexpi({imf.name}, 'L(T5|T4|E7|C8|ND)(\w*)', 'match');
imf = [imf{:}];
imf = vertcat(imf{:});
% sort according to yeardoy
yeardoy = str2num(imf(:, 10:16));
[~, sort_order] = sort(yeardoy);
imf = imf(sort_order, :);
% name of the first stacked image
filename = dir([imf(1,:),'/','L*MTLstack']); 
% read in dimension and zone number of the data
[jiDim,jiUL,resolu,zc] = envihdrread([imf(1,:),'/',filename.name]);
% dimension of image [row,col]
ijdim = [jiDim(2),jiDim(1)];
% number of nrows processed
nrows = ijdim(1);
% number of pixels procesed per line
ncols = ijdim(2);
% total ouput bands (1-5,7,6,cfmask)
nbands = 8;

rmse = zeros(nrows,ncols,nbands-1);

parfor i = 1:nrows
    i
    rmse(i,:,:) = Threshold_Line(ncols,i,nbands);
end

enviwrite('rmse',rmse,'int16',resolu,jiUL,'bsq',zc);

end % end of function 