
% load the saved tforms, the output format is good to copy and paste into
% my python script.
clc;
for iE = 1:13
    iB = iE + 1;
    try
        A = tforms{iB}.T';
    catch
        A = tforms{iB}.tdata.T';
    end

    M = zeros(4,4);
    M(1:2,1:2) = A(1:2,1:2);
    M(3,3) = 1;
    M(4,4) = 1;
    M(4,1:2) = A(3,1:2);
    M(1:2,4) = A(1:2,3);
    M;
    
    disp(['matrix: iE=',num2str(iE)]);
    
    str = sprintf('\t[%.4f, %.4f, 0, %.4f],\n\t[%.4f, %.4f, 0, %.4f],\n\t[0, 0, 1, 0],\n\t[0, 0, 0, 1]', ...
        M(1,1), M(1,2), M(1,4), M(2,1), M(2,2), M(2,4));
    
    disp(str);
    disp('===================================');
end