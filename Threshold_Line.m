function adj_rmse = Threshold_Line(ncols,nrows,nbands)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  defining variables
%% Constants
% number of bytes: int16
num_byte = 2;

% get num of total folders start with "L"
imf = dir('L*'); % folder names
% filter for Landsat folders
imf = regexpi({imf.name}, 'L(T5|T4|E7|C8|ND)(\w*)', 'match');
imf = [imf{:}];
imf = vertcat(imf{:});
% sort according to yeardoy
yeardoy = str2num(imf(:, 10:16));
[~, sort_order] = sort(yeardoy);
imf = imf(sort_order, :);
% number of folders start with "L"
num_t = size(imf,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Get ready for Xs & Ys
%% Read in Xs & Ysq
% transforming to serial date number (0000 year)
sdate = zeros(num_t,1); % Xs
line_t = zeros(num_t,nbands*ncols); %Ys
adj_rmse = zeros(ncols,nbands-1);
for i=1:num_t
    im_dir = dir(imf(i, :));
    im = '';
    for f = 1:size(im_dir, 1)
        % use regular expression to match:
        %   'L(\w*)'    Any word begining with L that has any following chars
        %   stk_n       includes stack name somewhere after L
        %   '$'           ends with the stack name (e.g., no .hdr, .aux.xml)
        if regexp(im_dir(f).name, ['L(\w*)', 'stack', '$']) == 1
            im = [imf(i, :), '/', im_dir(f).name];
            break
        end
    end
    % Check to make sure we found something
    if strcmp(im, '')
        error('Could not find stack image for directory %s\n', imf(i));
    end
    % Find date for folder imf(i)
    yr = str2num(imf(i, 10:13));
    doy = str2num(imf(i, 14:16));
    sdate(i) = datenum(yr, 1, 0) + doy;
    dummy_name = im;
    fid_t = fopen(dummy_name,'r'); % get file ids
    fseek(fid_t,num_byte*(nrows-1)*ncols*nbands,'bof');
    line_t(i,:) = fread(fid_t,nbands*ncols,'int16=>double','ieee-le'); % get Ys
end
fclose('all'); % close all files

for i_ids = 1:ncols
    % mask data
    line_m = line_t(:,nbands*i_ids);
    
    % Only run CCDC for places where more than 50% of images has data
    idexist = line_m < 255;
    overlap_pct = sum(idexist)/num_t;
    if overlap_pct < 0.5
        continue;
    end
    
    % convert Kelvin to Celsius
    line_t(:,nbands*(i_ids-1)+7) = line_t(:,nbands*(i_ids-1)+7)*10 - 27315;
    
    % clear pixel should have reflectance between 0 and 1
    % brightness temperature should between -93.2 to 70.7 celsius degree
    idrange = line_t(:,nbands*(i_ids-1)+1)>0&line_t(:,nbands*(i_ids-1)+1)<10000&...
        line_t(:,nbands*(i_ids-1)+2)>0&line_t(:,nbands*(i_ids-1)+2)<10000&...
        line_t(:,nbands*(i_ids-1)+3)>0&line_t(:,nbands*(i_ids-1)+3)<10000&...
        line_t(:,nbands*(i_ids-1)+4)>0&line_t(:,nbands*(i_ids-1)+4)<10000&...
        line_t(:,nbands*(i_ids-1)+5)>0&line_t(:,nbands*(i_ids-1)+5)<10000&...
        line_t(:,nbands*(i_ids-1)+6)>0&line_t(:,nbands*(i_ids-1)+6)<10000&...
        line_t(:,nbands*(i_ids-1)+7)>-9320&line_t(:,nbands*(i_ids-1)+7)<7070;
    
    % # of clear observatons
    idclr = line_m < 2;
    
    % clear and within physical range pixels
    idgood = idclr & idrange;
    
    % Xs & Ys for computation
    clrx = sdate(idgood);
    % bands 1-5,7,6
    clry = line_t(idgood,(nbands*(i_ids-1)+1):(nbands*(i_ids-1)+nbands-1));
    
    % caculate median variogram
    var_clry = clry(2:end,:)-clry(1:end-1,:);
    adj_rmse(i_ids,:) = median(abs(var_clry),1);
    
    
end % end of for i_ids loop
end % end of function
