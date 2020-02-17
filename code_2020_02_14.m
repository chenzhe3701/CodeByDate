% 2020-02-14
% Try RDR for a small sample area

% How to choose nDataPtsRange ========================================== 
% Assume DIC subsetSize = 2k+1, stepSize = s, filterSize = 2f+1
% A subset centered at ceil(k/s) will not be affected by the slipTraceLine 
% Due to averaging by filter, the non-affected subset center extended away from the slipTraceLine further by f data points
% Therefore, the non-affected subset center is ceil(k/s)+f data points away, or k+s*f pixels away. 
% In this study, nDataPtsRange > ceil(12/2)+2, so 10 seems fine.
nDataPtsRange = 10; % number of data points to cover on the horizontal line  

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
% Select area of interest (r4c5, r7c5, ...) ============================
<<<<<<< HEAD
ir=7;
=======
ir=4;
>>>>>>> origin/master
ic=5;
filename = ['E:\Ti7Al_E1_insitu_tension\Selected Area\r',num2str(ir),'c',num2str(ic),'\Ti7Al_E1_S7_r',num2str(ir),'c',num2str(ic),'.mat']
load(filename, 'x','y','exx','u','v','sigma');
exx(sigma==-1) = nan;
u(sigma==-1) = nan;
v(sigma==-1) = nan;

%% Illustrate ID map
myplot(X,Y,ID,boundaryTF);
label_map_with_ID(X,Y,ID,gcf,gID);

%% Plot map. For select grain (e.g., grain 9), predict theoretical RDR, and show selected slip traces 
<<<<<<< HEAD
ID_current = 226;
=======
ID_current = 176;
>>>>>>> origin/master
ind = (gID==ID_current);
euler = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
[ssa, c_a, nss, ntwin, ssGroup] = define_SS('Ti','pyii'); % [plane normal; slip direction] 
[ss, c_a, ssa] = define_SS_cart('Ti','pyii');
m = angle2dcm(euler(1)/180*pi, euler(2)/180*pi, euler(3)/180*pi, 'zxz');

clear traceDir bDir RDR;
for iss = 1:nss
    % trace direction for ref
    t = cross((ss(1,:,iss) * m), [0 0 1]);
    t = t(1,1:2);
    traceDir(iss,:) = t./norm(t);
    % RDR
    bDir(iss,:) = ss(2,:,iss) * m;
    rdr_t(iss,:) = bDir(iss,1)/bDir(iss,2); % theoretical, du/dv
end

% Plot map, draw theoretical slip traces for this grain to compare
<<<<<<< HEAD
ss_to_compare = [4,5];  % selected slip systems to draw trace and compare  
=======
ss_to_compare = [4,6,26];  % selected slip systems to draw trace and compare  
>>>>>>> origin/master
close all;

myplot(x,y,exx);
colormap(summer);
stressTensor = [1 0 0; 0 0 0; 0 0 0];
sampleMaterial = 'Ti';
label_map_with_trace_for_Ti7Al_E1(x,y,ones(size(x))*ID_current,ID_current,ss_to_compare,gca);    
tbl = array2table([[1:nss]',atand(traceDir(:,2)./traceDir(:,1)),rdr_t(:)]);
tbl.Properties.VariableNames = {'ss#','traceDir','rdr_theoretical'};
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
save(['E:\Ti7Al_E1_insitu_tension\Selected Area\selected_slip_trace_line_r',num2str(ir),'c',num2str(ic),'.mat'],'savedLines');
%% load savedLines
load(['E:\Ti7Al_E1_insitu_tension\Selected Area\selected_slip_trace_line_r',num2str(ir),'c',num2str(ic),'.mat'],'savedLines');
%% Load all saved lines for analysis, and repeat the above analysis to generate a map of rdr for each line on the slip trace line  

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
    text(pos(1,1)+100,pos(1,2),['Line ',num2str(iLine),', rdr=',num2str(rdr_exp,3)],'color','r');
    rdrLineMap(inds) = rdr_exp;
    rdrLines(iLine) = rdr_exp;
end
%% can plot the pt map and line map of rdr
myplot(x,y,rdrPtMap);
myplot(x,y,rdrLineMap);


%% repeat for each strain level, get rdr evolution at different strain levels  
clear rdrLine;
for iE = 1:7
    filename = ['E:\Ti7Al_E1_insitu_tension\Selected Area\r4c5\Ti7Al_E1_S',num2str(iE),'_r4c5.mat'];
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
<<<<<<< HEAD

        end
        lmd = fitlm(uv2(2,:),uv2(1,:));
        rdr_exp = lmd.Coefficients{2,1};
        
        rdrLine{iLine}(iE) = rdr_exp;
    end

end

=======

        end
        lmd = fitlm(uv2(2,:),uv2(1,:));
        rdr_exp = lmd.Coefficients{2,1};
        
        rdrLine{iLine}(iE) = rdr_exp;
    end

end

>>>>>>> origin/master
% plot
figure; hold on;
for iLine = 1:length(rdrLine)
    plot(rdrLine{iLine},'-o');
end
ylabel('rdr (experimental)');
xlabel('strain level');
legend('Location','bestoutside');




