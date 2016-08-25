function ARD_enviwrite(filename,envi_data,datatype,interleave,inputname)
%  write envi format file of uint8 and int16 
%  for example: enviwrite('testzz',data,'int16',[30,30],jiUL,'bsq',ZC,'this is for test only');
%  i,j => x, y => rows and cols
%  Input  1)  'filename' is the name of the file
%  Input  2)  'envi_data' is the data matrix
%  Input  3)  'dataype' is the datatype to write 'uint8' or 'int16'
%  Input  4)  'resolu' is the resolution of the pixel
%  Input  5)  'UL' is the UpperLeftPoint X Y of the UL pixel (not center of UL pixel)
%  Input  6)  'interleave' is the bsq bip bil format
%  Input  7)  'zone' is the UTM zone

if strcmp(datatype,'uint8')
    envi_data=uint8(envi_data);
    dt=1;
elseif strcmp(datatype,'int16')
    envi_data=int16(envi_data);
    dt=2;
elseif strcmp(datatype,'int32')
    envi_data=int32(envi_data);
    dt=3;
elseif strcmp(datatype,'uint16')
    envi_data=uint16(envi_data);
    dt=12;
elseif strcmp(datatype,'single')
    envi_data=single(envi_data);
    dt=4;
else
    fprintf('Invalid write data type!\n');
    return;
end

n_dims=size(envi_data);

nrows=n_dims(1); 
ncols=n_dims(2); 
bands=1;

if length(n_dims)>=3
    bands=n_dims(3);
end

multibandwrite(envi_data,filename,interleave);

filename_hdr=[filename,'.hdr'];
                                
fid_out=fopen(filename_hdr,'wt');

fprintf(fid_out,'ENVI\n');
fprintf(fid_out,'descirption = {Landsat Scientific Data}\n');
fprintf(fid_out,'samples = %d\n',ncols); % samples is for j
fprintf(fid_out,'lines   = %d\n',nrows); % lines is for i
fprintf(fid_out,'bands   = %d\n',bands);
fprintf(fid_out,'header offset = 0\n');
fprintf(fid_out,'file type = ENVI Standard\n');
fprintf(fid_out,'data type = %d\n',dt);
fprintf(fid_out,'interleave = %s\n',interleave);
fprintf(fid_out,'sensor type = Landsat\n');
fprintf(fid_out,'byte order = 0\n');

% Read metadata from example_img.hdr or HDR
filename_HDR=[inputname,'.HDR'];
filename_hdr=[inputname,'.hdr'];

fid_in1=fopen(filename_hdr,'r');
fid_in2=fopen(filename_HDR,'r');

if fid_in1~=-1
    fid_in=fid_in1;
elseif fid_in2~=-1
    fid_in=fid_in2;
else
    fprintf('Wrong input ENVI header file!\n');
    sprintf('%s\n',filename); % show envi hdr file
    return;
end

geo_char=fscanf(fid_in,'%c',inf);

map_pos = strfind(geo_char,'map info = ');
proj_pos = strfind(geo_char,'projection info = ');
coor_pos = strfind(geo_char,'coordinate system string = ');
band_pos = strfind(geo_char,'band names = ');

fprintf(fid_out, '%s',geo_char(map_pos:proj_pos-1));
fprintf(fid_out, '%s',geo_char(proj_pos:coor_pos-1));
fprintf(fid_out, '%s',geo_char(coor_pos:band_pos-1));

fclose(fid_out);
end



