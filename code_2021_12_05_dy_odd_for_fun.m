
n_to_success = [];
prob = 0.119;

N = round(10000/12);    % total number of trials

ii = 1;
count = 1;

while ii < N
    n = rand(1);
    if n<prob
        r = true;
        n_to_success(ii) = count;
        count = 1;
    else
        r = false;
        count = count + 1;
    end
    
    ii = ii + 1;
end

n_to_success(n_to_success==0) = [];
max(n_to_success)
s = sort(n_to_success, 'descend');