function ShowLCMapCrop
% This function is used to provde all change maps for each year
% Results for LCMAP

% INPUTS:

v_input = ccdc_Inputs;
pwd

% read dem as mask
dem = enviread([pwd,'/ANC/dem']);

% load LCMap products
ChangeMap = enviread([pwd,'/FullCCDCMap/ChangeMap']);
ChangeMagMap = enviread([pwd,'/FullCCDCMap/ChangeMagMap']);
CoverMap = enviread([pwd,'/FullCCDCMap/CoverMap']);
CoverQAMap = enviread([pwd,'/FullCCDCMap/CoverQAMap']);
ConditionMap = enviread([pwd,'/FullCCDCMap/ConditionMap']);
QAMap = enviread([pwd,'/FullCCDCMap/QAMap']);

% dimension and projection of the image
nrows = v_input.ijdim(1);
ncols = v_input.ijdim(2);
jiDim = [ncols,nrows];
jiUL = v_input.jiul;
res = v_input.resolu;
zc = v_input.zc;
nbands = size(QAMap,3);
% make Predict folder for storing predict images
n_map = v_input.name_map;% 'CCDCMap';

for i = 1:nbands
    
    i
    
    TMP = ChangeMap(:,:,i);
    TMP(dem == 0) = 0;
    ChangeMap(:,:,i) = TMP;
    
    TMP = ChangeMagMap(:,:,i);
    TMP(dem == 0) = 0;
    ChangeMagMap(:,:,i) = TMP;
    
    TMP = CoverMap(:,:,i);
    TMP(dem == 0) = 0;
    CoverMap(:,:,i) = TMP;
    
    TMP = CoverQAMap(:,:,i);
    TMP(dem == 0) = 0;
    CoverQAMap(:,:,i) = TMP;
    
    TMP = ConditionMap(:,:,i);
    TMP(dem == 0) = 0;
    ConditionMap(:,:,i) = TMP;
    
    TMP = QAMap(:,:,i);
    TMP(dem == 0) = 0;
    QAMap(:,:,i) = TMP;
end  

% write cropped LCMap products
% write ENVI files
% Cover Map
enviwrite([v_input.l_dir,'/',n_map,'/CoverMap'],CoverMap,'uint8',res,jiUL,'bsq',zc);
enviwrite([v_input.l_dir,'/',n_map,'/CoverQAMap'],CoverQAMap,'uint8',res,jiUL,'bsq',zc);

% Change Map
enviwrite([v_input.l_dir,'/',n_map,'/ChangeMap'],ChangeMap,'uint16',res,jiUL,'bsq',zc);
enviwrite([v_input.l_dir,'/',n_map,'/ChangeMagMap'],ChangeMagMap,'single',res,jiUL,'bsq',zc);

% Condition Map
enviwrite([v_input.l_dir,'/',n_map,'/ConditionMap'],ConditionMap,'single',res,jiUL,'bsq',zc);
% Model Estimation QA Map
enviwrite([v_input.l_dir,'/',n_map,'/QAMap'],QAMap,'uint8',res,jiUL,'bsq',zc);


