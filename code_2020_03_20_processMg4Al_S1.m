% scripts used to process Mg4Al_S1 sample

addChenFunction;
cd('E:\Mg4Al_S1_insitu\SEM_Images\stitched_DIC');

%% plot strain map
for ii=1:6
   load(['global_method_stitched_',num2str(ii),'.mat'],'exx');
   myplot(exx);
   print(['exx_iE_',num2str(ii),'.tif'],'-dtiff');
end

%% plot a map to identify control points
ii = 3
load(['global_method_stitched_',num2str(ii),'.mat'],'x','y','exx');
myplot(x,y,exx);

cpEBSD = [250, 280;
    1247, 268;
    280, 1186;
    1283, 1307;
    
    288, 218;
    1231, 156;
    357, 1415;
    1239, 1425] / 2;    % pixel location/2 -> um location

cpSEM = [1820, 5520;
    41980, 4800;
    3540, 37620;
    43380, 41060;
    
    3340, 3400;
    41400, 780;
    6740, 45580;
    41660, 45200];