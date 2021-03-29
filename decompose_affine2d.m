% Decompose affine transformation using python library transforms3d
% varInput can be 
% (1) the tform from maketform('affine',cpFrom,cpTo)
% (2) the tform from fitgeotrans(cpFrom, cpTo, 'affine')
% (3) or the matrix A for transformation, where cpTo_3x1 = A * cpFrom_3x1
% output: T=translation, R=rotation matrix, Z=zoom, S=shear
% chenzhe, 2021-03-10

function [T,R,Z,S] = decompose_affine2d(varInput)

% add library path
setenv('path',['C:\Users\ZheChen\miniconda3\envs\pythonForMatlab\Library\bin;', getenv('path')]); 
pe = pyenv;
if strcmpi(pe.Status,'NotLoaded')
    setenv('PYTHONUNBUFFERED','1'); % doesn't change much, but just-in-case 
    pe = pyenv('Version','C:\Users\ZheChen\miniconda3\envs\pythonForMatlab\python.exe',"ExecutionMode","InProcess");
    % py.importlib.import_module('numpy')
end

if isnumeric(varInput)
    A2 = varInput;
elseif isobject(varInput)
    A2 = varInput.T';
elseif isstruct(varInput)
    A2 = varInput.tdata.T';
else
    error('unknown data type');
end

A = [A2(1,1), A2(1,2), 0, A2(1,3);
     A2(2,1), A2(2,2),  0, A2(2,3);
     0,       0,        1, 0;
     A2(3,1), A2(3,2),  0, A2(3,3)];

result = py.transforms3d.affines.decompose(A);   % result = T,R,Z,S

T = double(result{1});
R = double(result{2});
Z = double(result{3});
S = double(result{4});

end