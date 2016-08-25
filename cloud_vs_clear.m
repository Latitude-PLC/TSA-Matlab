% Start with Boston images
% pct_mat1=test('/net/casfsd/fs/scratch11/zhuzhe/');
% pct_mat2=test('/net/casfsd/fs/scratch12/zhuzhe/');
% pct_t=test('/net/casfsd/fs/scratch13/zhuzhe/');
% pct_t=test('/net/casfsd/fs/scratch31/Lp17r37/');
% pct_t=[pct_mat1,pct_mat2];
% start histograming

figure;
nsize = 14;
set(gca,'FontSize',nsize);
x = 0:0.025:1;
nbar = length(x);
y = zeros(1,nbar);
num = zeros(1,nbar);
total_clr = sum(prct)/100;
pct_t = prct/100;

for i=1:nbar-1
    l_pct = length(pct_t(pct_t<x(i+1)&pct_t>=x(i)));
    num(i) = l_pct;
    y(i) = sum(pct_t(pct_t<x(i+1)&pct_t>=x(i)));
    if i>1
        y(i) = y(i-1)+y(i);
        num(i) = num(i-1)+num(i);
    end
end


x=100*x;
hist_prct=100*y/y(nbar-1);
hist_num = 100*num/num(nbar-1);

% bar(x(1:end-1),hist_prct(1:end-1),'b');
% hold on
plot(x(1:end-1),hist_num(1:end-1),'r');
hold on


bar(x,hist_prct,'hist');
% for i=1:nbar-1
%     text(x(i),103,num2str(num(i)),'FontSize',nsize);
% end
%bar(100*x,num,'r','hist');
xlabel('Clear Percent Interval (%)');
ylabel('Percent of Total Clear Observations (%)');
axis([0 100 0 100])
% title({'Massachusetts, USA';''});
%title({'Georgia, USA';''});
 title({'grid02';''});
grid on;

% % function pct_mat=test(dir_l)
% % % locations of data
% % % cd /net/casfsd/fs/scratch14/zhuzhe/ % Boston scene
% % cd(dir_l);
% % % Variables
% % % % get num of total folders start with "L"
% % imf=dir('L*'); % folder names
% % num_t=size(imf,1);
% % obs_n='fmask_1';
% % pct_mat=zeros(1,num_t);
% % 
% % for i=1:num_t
% %    fprintf('Processing the %dth image\n',i);
% %    im=auto_imget(i,'L');% get image by order of the folder name
% %    Fmask=enviread([im.n_name,obs_n]); % read in whole image
% %    cld_pix=sum(sum(Fmask==3)); % num of cloud pixels
% %    others_pix=sum(sum(Fmask<=3)); % number of others pixels (including clouds)
% %    pct_mat(i)=cld_pix/others_pix;
% % end
