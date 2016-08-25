function MainClassMap(yy,mm,dd)
% This function is used to provde Change Type Maps for each pixels

% yy is a matrix of year
% mm is the month
% dd is the day

for i = 1:length(yy);
    fprintf('Classification map at %d/%d/%d is processing ...!\n',mm,dd,yy(i));
    ShowClassMap(yy(i),mm,dd);
end