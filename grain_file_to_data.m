% 2021-01-24
% load from f1 = grain file type 1, f2 = grain file type 2
function data = grain_file_to_data(f1,f2)

try
    [EBSD_data_1, EBSD_header_1] = grain_file_read(f1);
    % find column
    column_index_1 = find_variable_column_from_grain_file_header(EBSD_header_1, ...
        {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
    
    % EBSD data, from type-1 grain file. (column, data) pair:
    % (1,phi1) (2,phi) (3,phi2) (4,xMicron) (5,yMicron) (6,IQ) (7,CI) (8,Fit) (9,grain-Id) (10,edgeGrain?)
    % Read EBSD data.  IQ,CI,Fit are not needed for now, but might need in future
    x = EBSD_data_1(:,column_index_1(5));
    y = EBSD_data_1(:,column_index_1(6));
    unique_x = unique(x(:));
    ebsdStepSize = unique_x(2) - unique_x(1);
    mResize = (max(x(:)) - min(x(:)))/ebsdStepSize + 1;
    nResize = (max(y(:)) - min(y(:)))/ebsdStepSize + 1;
    
    data.x = reshape(EBSD_data_1(:,column_index_1(5)),mResize,nResize)';
    data.y = reshape(EBSD_data_1(:,column_index_1(6)),mResize,nResize)';
    data.ID = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
    phi1 = reshape(EBSD_data_1(:,column_index_1(2)),mResize,nResize)';
    phi = reshape(EBSD_data_1(:,column_index_1(3)),mResize,nResize)';
    phi2 = reshape(EBSD_data_1(:,column_index_1(4)),mResize,nResize)';
    % change it to degrees, if necessary
    if max(phi1(:))<7 && max(phi(:))<7 && max(phi2(:))<7
        phi1 = phi1*180/pi();
        phi = phi*180/pi();
        phi2 = phi2* 180/pi();
    end
    data.phi1 = phi1;
    data.phi = phi;
    data.phi2 = phi2;
catch
end

try
    [EBSD_data_2, EBSD_header_2] = grain_file_read(f2);
    column_index_2 = find_variable_column_from_grain_file_header(EBSD_header_2, ...
        {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});
    
    % read type-2 grain file and get average info for grains
    data.gID = EBSD_data_2(:,column_index_2(1));
    data.gPhi1 = EBSD_data_2(:,column_index_2(2));
    data.gPhi = EBSD_data_2(:,column_index_2(3));
    data.gPhi2 = EBSD_data_2(:,column_index_2(4));
    data.gCenterX = EBSD_data_2(:,column_index_2(5));
    data.gCenterY = EBSD_data_2(:,column_index_2(6));
    data.gEdge = EBSD_data_2(:,column_index_2(10));
    
    data.gNNeighbors = EBSD_data_2(:,column_index_2(7));
    data.gDiameter = EBSD_data_2(:,column_index_2(8));
    data.gArea = EBSD_data_2(:,column_index_2(9));
    data.gNeighbors = EBSD_data_2(:,(column_index_2(7)+1):(size(EBSD_data_2,2)));
catch
end