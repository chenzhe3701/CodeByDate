% 2020-02-14
% Try RDR for a small sample area

% Grain #1-9 in Hank's map have:
% ID = [80, 226, 185, 221, 83, 95, 189, 130, 176] in my EBSD data
clear;clc;
addChenFunction;
addpath('D:\p\m\CodeByDate');

% How to choose nDataPtsRange ========================================== 
% Assume DIC subsetSize = 2k+1, stepSize = s, filterSize = 2f+1
% A subset centered at ceil(k/s) will not be affected by the slipTraceLine 
% Due to averaging by filter, the non-affected subset center extended away from the slipTraceLine further by f data points
% Therefore, the non-affected subset center is ceil(k/s)+f data points away, or k+s*f pixels away. 
% In this study, nDataPtsRange > ceil(12/2)+2, so 10 seems fine.
nDataPtsRange = 10; % number of data points to cover on the horizontal line  

% Load orientation data
load('E:\Ti7Al_E1_insitu_tension\Analysis_by_Matlab\Ti7Al_E1_organized.mat', 'gID','gPhi1','gPhi','gPhi2','ID','eulerAligned'); 
load('E:\Ti7Al_E1_insitu_tension\Analysis_by_Matlab\Ti7Al_E1_EbsdToSemForTraceAnalysis.mat', 'X','Y','boundaryTF')
load('D:\p\m\DIC_Analysis\setting_for_real_samples\Ti7Al_E1_setting.mat','strainPauses_DIC');
% Select area of interest (r4c5, r7c5, r6c5 ...) ============================
ir = 6;
ic = 5;
iE = 7;
filename = ['E:\Ti7Al_E1_insitu_tension\Analysis_selected_area\r',num2str(ir),'c',num2str(ic),'\Ti7Al_E1_S',num2str(iE),'_r',num2str(ir),'c',num2str(ic),'.mat']
load(filename, 'x','y','exx','u','v','sigma');
exx(sigma==-1) = nan;
u(sigma==-1) = nan;
v(sigma==-1) = nan;

exx_all = load('E:\Ti7Al_E1_insitu_tension\SEM Data\stitched_DIC\global_method_stitched_7.mat','exx');
exx_all = exx_all.exx;
xum = X / 4096 * 60;
yum = Y / 4096 * 60;

%% Illustrate ID map
myplot(X,Y,ID,boundaryTF);
label_map_with_ID(X,Y,ID,gcf,gID);

%% illustrate exx map
[f,a,c] = myplot(xum,yum,exx_all,boundaryTF);
title('');
set(gca,'XTick', [0:100:400], 'YTick', [0:100:400], 'fontsize',18);
xlabel('x (\mum)');
ylabel('y (\mum)');
title(c,'\epsilon_x_x');

%% Can load local boundary map
load('E:\Ti7Al_E1_insitu_tension\Analysis_r6c5\Ti7Al_E1_EbsdToSemForTraceAnalysis.mat','boundaryTF')

%% Plot map. For select grain (e.g., grain 9), predict theoretical RDR, and show selected slip traces 
ID_current = 226;
ind = (gID==ID_current);
euler = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
[ssa, c_a, nss, ntwin, ssGroup] = define_SS('Ti','pyii'); % [plane normal; slip direction] 
[ss, c_a, ssa] = define_SS_cart('Ti','pyii');
m = angle2dcm(euler(1)/180*pi, euler(2)/180*pi, euler(3)/180*pi, 'zxz');
stressTensor = [1 0 0; 0 0 0; 0 0 0];

clear traceDir bDir RDR;
for iss = 1:nss
    % trace direction for ref
    t = cross((ss(1,:,iss) * m), [0 0 1]);
    t = t(1,1:2);
    traceDir(iss,:) = t./norm(t);
    % RDR
    bDir(iss,:) = ss(2,:,iss) * m;
    rdr_t(iss,:) = bDir(iss,1)/bDir(iss,2); % theoretical, du/dv
    
    N = ss(1,:,iss) * m;    % plane normal in sample coord
    M = ss(2,:,iss) * m;    % slip dir in sample coord
    sf(iss,:) = abs(N * stressTensor * M');
end

% Plot map, draw theoretical slip traces for this grain to compare
ss_to_compare = [4,5];  % selected slip systems to draw trace and compare  
close all;

myplot(x,y,exx,boundaryTF);
colormap(summer);
stressTensor = [1 0 0; 0 0 0; 0 0 0];
sampleMaterial = 'Ti';
label_map_with_trace_for_Ti7Al_E1(x,y,ones(size(x))*ID_current,ID_current,ss_to_compare,'pyii',gca);    
tbl = array2table([[1:nss]', atand(traceDir(:,2)./traceDir(:,1)), rdr_t(:), sf(:)]);
tbl.Properties.VariableNames = {'ss#','traceDir','rdr_theoretical','sf'};
disp(tbl);

%% Add line on top of slip trace to analyze
handleDrawline = drawline(gca,'Color','r');
pos = customWait(handleDrawline)

% Analyze data
% find coordinate/maybe indices
pt1 = pos(1,:);
pt2 = pos(2,:);

% find intersection between the line drawn and uniqueGrainBoundary.
[pts,inds,indR,indC] = grids_covered_by_line(x,y,pt1,pt2);

uv = [];
uv2 = [];
uvAdjusted = [];
uv_range = [];
tAngle = atand((pt2(2)-pt1(2))/(pt2(1)-pt1(1)));   % Angle of the labeled slip trace line. assume line is not verticle.
pAngle = tAngle + 90;   % angle of the perpendicular line of slip trace line.
indRange = ceil(max(abs([sind(pAngle),cosd(pAngle)]))*nDataPtsRange);

figure; hold on;
for ii = 1:length(inds)
    indr = indR(ii);    % index of the current data point on slip trace line
    indc = indC(ii);
    
    xLocal = x(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
    yLocal = y(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
    uLocal = u(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
    vLocal = v(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
    exxLocal = exx(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
    sigmaLocal = sigma(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
    
    [ptLocal,indLocal,indrLocal,indcLocal] = grids_covered_by_line(xLocal, yLocal, pts(ii,:)-[1000,1000*tand(pAngle)], pts(ii,:)+[1000,1000*tand(pAngle)]);

    uvLocal = [uLocal(indLocal),vLocal(indLocal)]'; % row 1 u, row 2 v
    % old-method, consider all data points on perpendicular line
    %     row_mean = nanmean(uvLocal,2);
    %     uvLocal = uvLocal - repmat(row_mean,1,size(uvLocal,2));      % displacements u,v, and post-deformed locations x_p,y_p are all centered.
    %     uv = [uv,uvLocal];
    %     uv_range = [uv_range,max(uvLocal,[],2)-min(uvLocal,[],2)];
    %     plot(uvLocal(1,:),uvLocal(2,:),'-o','color', rand(1,3));

    % new method, consider only two end points on perpendicular line
    uvLocal = uvLocal(:,[1,end]);
    row_mean = nanmean(uvLocal,2);
    uvLocal = uvLocal - repmat(row_mean,1,size(uvLocal,2)); 
    uv2 = [uv2,uvLocal];
    plot(uvLocal(1,:),uvLocal(2,:),'-dr','linewidth',2)

end
lmd = fitlm(uv2(2,:),uv2(1,:));
rdr_exp = lmd.Coefficients{2,1}

%% Can append this slip trace line to saved slip trace lines for analysis
if ~exist('savedLines','var')
    savedLines = {};
end
nL = length(savedLines);
savedLines{nL+1} = pos;

%% Can plot u,v,sigma fields for illustration
myplot(x,y,u);
colormap(summer);
drawline(gca,'Position',pos);

myplot(x,y,v);
colormap(summer);
drawline(gca,'Position',pos);

myplot(x,y,sigma);
colormap(summer);
drawline(gca,'Position',pos);

%% save savedLines
save(['E:\Ti7Al_E1_insitu_tension\Analysis_selected_area\selected_slip_trace_line_r',num2str(ir),'c',num2str(ic),'.mat'],'savedLines');
%% load savedLines
load(['E:\Ti7Al_E1_insitu_tension\Analysis_selected_area\selected_slip_trace_line_r',num2str(ir),'c',num2str(ic),'.mat'],'savedLines');
%% Load all saved lines for analysis, and repeat the above analysis to generate a map of rdr for each line on the slip trace line  

myplot(x,y,exx,boundaryTF);title('\epsilon_x_x','fontweight','normal');
colormap(summer);
set(gca,'fontsize',18,'XTick',[0:1000:4000],'YTick',[0:1000:4000]);
xlabel('X (pixels)');
ylabel('Y (pixels)');

% This is for r6c5
ir_aoi = 1060;
ic_aoi = 280;
rectangle('position',[2*ic_aoi,2*ir_aoi,2*340,2*340],'linewidth',2);
label_map_with_trace_for_Ti7Al_E1(x,y,ones(size(x))*ID_current,ID_current,ss_to_compare,'pyii',gca);   


rdrLineMap = zeros(size(x));
rdrPtMap = zeros(size(x));
for iLine = 1:length(savedLines)
    pos = savedLines{iLine};
    drawline(gca,'Position',pos,'color','r')
    % find coordinate/maybe indices
    pt1 = pos(1,:);
    pt2 = pos(2,:);
    
    % find intersection between the line drawn and uniqueGrainBoundary.
    [pts,inds,indR,indC] = grids_covered_by_line(x,y,pt1,pt2);
    
    uv = [];
    uv2 = [];
    uvAdjusted = [];
    uv_range = [];
    tAngle = atand((pt2(2)-pt1(2))/(pt2(1)-pt1(1)));   % Angle of the labeled slip trace line. assume line is not verticle.
    pAngle = tAngle + 90;   % angle of the perpendicular line of slip trace line.
    indRange = ceil(max(abs([sind(pAngle),cosd(pAngle)]))*nDataPtsRange);
    
    for ii = 1:length(inds)
        indr = indR(ii);    % index of the current data point on slip trace line
        indc = indC(ii);
        
        xLocal = x(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
        yLocal = y(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
        uLocal = u(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
        vLocal = v(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
        exxLocal = exx(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
        sigmaLocal = sigma(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
        
        [ptLocal,indLocal,indrLocal,indcLocal] = grids_covered_by_line(xLocal, yLocal, pts(ii,:)-[1000,1000*tand(pAngle)], pts(ii,:)+[1000,1000*tand(pAngle)]);
        
        uvLocal = [uLocal(indLocal),vLocal(indLocal)]'; % row 1 u, row 2 v
        
        % new method, consider only two end points on perpendicular line
        uvLocal = uvLocal(:,[1,end]);
        row_mean = nanmean(uvLocal,2);
        uvLocal = uvLocal - row_mean;
        uv2 = [uv2,uvLocal];
        rdrPtMap(inds(ii)) = (uvLocal(1,2)-uvLocal(1,1))/(uvLocal(2,2)-uvLocal(2,1));    % (u2-u1)/(v2-v1) for perpendicular line, assigned to center data point
    end
    lmd = fitlm(uv2(2,:),uv2(1,:));
    rdr_exp = lmd.Coefficients{2,1};
%     text(pos(1,1)+100,pos(1,2),['Line ',num2str(iLine),', rdr=',num2str(rdr_exp,3)],'color','r','fontsize',18);
    text(pos(1,1)+100,pos(1,2),['Line ',num2str(iLine)],'color','r','fontsize',18);
    rdrLineMap(inds) = rdr_exp;
    rdrLines(iLine) = rdr_exp;
end
text(500,500,'Grain #1','color','k','fontsize',18);
%% can plot the pt map and line map of rdr
myplot(x,y,rdrPtMap);
myplot(x,y,rdrLineMap);


%% repeat for each strain level, get rdr evolution at different strain levels  
clear rdrLine;
for iE = 1:7
    filename = ['E:\Ti7Al_E1_insitu_tension\Analysis_selected_area\r',num2str(ir),'c',num2str(ic),'\Ti7Al_E1_S',num2str(iE),'_r',num2str(ir),'c',num2str(ic),'.mat']
    load(filename, 'x','y','exx','u','v','sigma')
    exx(sigma==-1) = nan;
    u(sigma==-1) = nan;
    v(sigma==-1) = nan;
    
    for iLine = 1:length(savedLines)
        pos = savedLines{iLine};
        % find coordinate/maybe indices
        pt1 = pos(1,:);
        pt2 = pos(2,:);
        
        % find intersection between the line drawn and uniqueGrainBoundary.
        [pts,inds,indR,indC] = grids_covered_by_line(x,y,pt1,pt2);
        
        uv = [];
        uv2 = [];
        uvAdjusted = [];
        uv_range = [];
        tAngle = atand((pt2(2)-pt1(2))/(pt2(1)-pt1(1)));   % Angle of the labeled slip trace line. assume line is not verticle.
        pAngle = tAngle + 90;   % angle of the perpendicular line of slip trace line.
        indRange = ceil(max(abs([sind(pAngle),cosd(pAngle)]))*nDataPtsRange);
        
        for ii = 1:length(inds)
            indr = indR(ii);    % index of the current data point on slip trace line
            indc = indC(ii);
            
            xLocal = x(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
            yLocal = y(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
            uLocal = u(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
            vLocal = v(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
            exxLocal = exx(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
            sigmaLocal = sigma(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
            
            [ptLocal,indLocal,indrLocal,indcLocal] = grids_covered_by_line(xLocal, yLocal, pts(ii,:)-[1000,1000*tand(pAngle)], pts(ii,:)+[1000,1000*tand(pAngle)]);
            
            uvLocal = [uLocal(indLocal),vLocal(indLocal)]'; % row 1 u, row 2 v
            
            % new method, consider only two end points on perpendicular line
            uvLocal = uvLocal(:,[1,end]);
            row_mean = nanmean(uvLocal,2);
            uvLocal = uvLocal - row_mean;
            uv2 = [uv2,uvLocal];

        end
        lmd = fitlm(uv2(2,:),uv2(1,:));
        rdr_exp = lmd.Coefficients{2,1};
        
        rdrLine{iLine}(iE) = rdr_exp;
        ebar{iLine}(iE) = std(uv2(1,:)./uv2(2,:));
%         ebar{iLine}(1:2) = nan;
    end

end

%% plot for r6c5
close all;figure; hold on;

% plot(1:7,rdrLine{1},'-o','linewidth',2);
% plot(1:7,rdrLine{2},'-^','linewidth',2);
% plot(1:7,rdrLine{3},'-d','linewidth',2);
% plot(1:7,rdrLine{4},'-s','linewidth',2);

errorbar(1:7,rdrLine{1},ebar{1},'-o','linewidth',2,'CapSize',12);
errorbar(1:7,rdrLine{2},ebar{2},'-^','linewidth',2,'CapSize',12);
errorbar(1:7,rdrLine{3},ebar{3},'-d','linewidth',2,'CapSize',12);
errorbar(1:7,rdrLine{4},ebar{4},'-v','linewidth',2,'CapSize',12);

plot(0:8, repmat(rdr_t(4),1,9), '--k','linewidth',2,'HandleVisibility','off')
text(3.5,rdr_t(4)-0.3,['RDR^4_{theoretical} = ',num2str(rdr_t(4),3)],'fontsize',16)
plot(0:8, repmat(rdr_t(5),1,9), '--k','linewidth',2,'HandleVisibility','off')
text(3.5,rdr_t(5)+0.3,['RDR^5_{theoretical} = ',num2str(rdr_t(5),3)],'fontsize',16)

set(gca,'ylim',[-2.2,2],'fontsize',18)
ylabel('RDR');
xlabel('Strain Level');
legend({'Line 1','Line 2','Line 3', 'Line 4', 'Line 5'},'Location','southeast');


%% plot for any new measurement
figure; hold on;
for iLine = 1:length(rdrLine)
    plot(1:length(rdrLine{iLine}),rdrLine{iLine},'-o','linewidth',2);
end
set(gca,'fontsize',18)
ylabel('RDR');
xlabel('Strain Level');
legend({'Line 1','Line 2','Line 3', 'Line 4', 'Line 5'},'Location','southeast');