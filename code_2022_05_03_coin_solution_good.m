% coin/ball problem, find a defect ball (do not know heavier or lighter) in
% a total of N balls, using a balance.
% The maximum N that can be processed in k measure is:
% N = (3^k-3)/2

for i = 1:1000
    clc;
    % make a vector of length N, value 1, randomly change one element
    N = randi([40,120]);

    v = ones(1, N);
    v(randi([1,N])) = rand + 0.5;
    

    defect_type = nan;    % -1 = smaller defect, 1 = larger defect, 0 = left portion smaller, but do not know relationship with good value.  
    count = 0;

    [defect, defect_type, count] = find_in_vector(v, defect_type, count);

    assert(defect ~= 1);
    assert(defect_type > -2);
    assert((defect-1) * defect_type > 0);
    assert(count <= 5);

    display(v);
    disp(['defect: ', num2str(defect)]);
    disp(['defect type: ', num2str(defect_type)]);
    disp(['counts: ', num2str(count)]);

end


% find defect in a vector
function [defect, defect_type, count] = find_in_vector(v, defect_type, count)

good_value = 1;
N = length(v);

if N == 1    
    defect = v(1);
    if isnan(defect_type)
        [status, count] = compare_vector(v(1), good_value, count);
        defect_type = status;
    end
elseif N == 2
    if abs(defect_type) == 1
        [status, count] = compare_vector(v(1),v(2), count);
        if status == defect_type
            defect = v(1);
        else
            defect = v(2);
        end
    elseif defect_type == 0
        [status, count] = compare_vector(v(1), good_value, count);
        if status == -1
            defect = v(1);
            defect_type = -1;
        elseif status == 0
            defect = v(2);
            defect_type = 1;
        elseif status == 1
            throw(MException('myfunction:error', 'logical error'));
        end
    elseif isnan(defect_type)
        [status, count] = compare_vector(v(1), good_value, count);
        if status == 0
            defect = v(2);
            [defect_type, count] = compare_vector(v(2),good_value,count);
        else
            defect = v(1);
            defect_type = status;
        end
    end
else
    if isnan(defect_type)
        [a,b,c,d,use_two_segments] = get_ab(N, defect_type, count);
        
        if use_two_segments == false        
            v1 = v(1:a);
            v2 = v(a+1 : 2*a);
            v3 = v(2*a+1 : end);
            [status, count] = compare_vector(v1, v2, count);
            
            if status == 0
                defect_type = defect_type;  % no updates on this
                v = v3;
            else
                if status == 1
                    defect_type = 0;
                    v = [v2,v1];
                else
                    defect_type = 0;
                    v = [v1,v2];
                end
            end
            [defect, defect_type, count] = find_in_vector(v, defect_type, count);
        
        elseif use_two_segments == true
            v1 = v(1:c);
            v2 = v(c+1:end);
            [status, count] = compare_vector(v1, ones(size(v1))*good_value, count);

            if status == 0
                v = v2;
            else
                defect_type = status;
                v = v1;
            end
            [defect, defect_type, count] = find_in_vector(v, defect_type, count);
        end
  
    elseif defect_type == 0
        % else, we know the defect type partially  
        % we divide into [A1 B1 | B2 A2], replace A1 with A2, A2 with good_values, and compare
        [a,b,c,d,use_two_segments] = get_ab(N, defect_type, count);
        assert(use_two_segments==false);
        
        [A1,v] = split_vector(v, a);
        [B1,v] = split_vector(v, b/2);
        [B2,v] = split_vector(v, b/2);
        [A2,v] = split_vector(v, a);
        assert(isempty(v));
        assert(rem(b,2)==0);    % if we partially know, b should always be even 

        [status, count] = compare_vector([A2,B1], [B2,ones(size(A2))*good_value], count);
        if status == -1
            v = [B1,B2];
        elseif status == 0
            % A2 is good, A1 is ligher
            defect_type = -1;
            v = A1;
        elseif status == 1
            % A2 is larger
            defect_type = 1;
            v = A2;
        end
        [defect, defect_type, count] = find_in_vector(v, defect_type, count);
        
    elseif abs(defect_type) == 1
        % else, we know there is only one possibility of the defect type  
        % we divide into [A1 B A2], compare A1 with A2
        [a,b,c,d,use_two_segments] = get_ab(N, defect_type, count);
        assert(use_two_segments==false);

        [A1,v] = split_vector(v,a);
        [B,v] = split_vector(v,b);
        [A2,v] = split_vector(v,a);
        assert(isempty(v));

        [status, count] = compare_vector(A1, A2, count);
        if status == 0
            % B has defect
            v = B;
        elseif status == defect_type
            % if order does not change, A1 has defect
            v = A1;
        else
            % if order changes, A2 has defect
            v = A2;
        end
        [defect, defect_type, count] = find_in_vector(v, defect_type, count);
    end
end

end


% compare 2 vectors v1 and v2, return -1, 0, or 1
function [status, count] = compare_vector(v1, v2, count)
if sum(v1) > sum(v2)
    status = 1;
elseif sum(v1) == sum(v2)
    status = 0;
else
    status = -1;
end
count = count + 1;
end


% split vector at index
function [va,vb] = split_vector(v,index)
va = v(1:index);
v(1:index) = [];
vb = v;
end


% get segment length
function [a,b,c,d,use_two_segments] = get_ab(N, defect_type, count)

a = nan;
b = nan;
c = nan;
d = nan;
use_two_segments = false;

if isnan(defect_type)
    % If we do not know what type of defect it is (e.g., initial setup [A1 A2 B], compare A1 A2, return [A1 A2] or [B])
    % N=2a+b, 2a<=2(b+2), 2b<=2(a+1) ==> a in [(N-1)/3, (N+2)/3], b in [(N-4)/3, (N+2)/3]
    a = ceil((N-1)/3);
    b = N - 2*a;
    
    % If count is larger than 0, that means we are gauranteed to have some good elements  
    % And if we still do not know anything about the defect type
    % Then it is possible that two segments can give a better solution...  
    if count > 0
        % if we already have some good element for reference, then we want to segment into [c,d], and compare c with good_values   
        % N=c+d, c+1>=2d, c<=2(d+1) ==> c in [(2N-1)/3, (2N+2)/3], d in [[(N-2)/3, (N+1)/3]
        c = ceil((2*N-1)/3);
        d = N - c;

        if max(2*a,2*b) > max(c, 2*d)
            use_two_segments = true;
        end
    end
else
    % If we know there is only one possible solution left for the defect type, 
    %  then we divide into [A1 B A2], replace A1 with A2, A2 with good_values, and compare    
    % N=2a+b, a<=b+2, b<=a+1 ==> a in [(N-1)/3, (N+2)/3], b in [(N-4)/3, (N+2)/3]   
    a = ceil((N-1)/3);
    b = N - 2*a;
end

end