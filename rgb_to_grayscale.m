% the input is rgb matrix
% the output is grayscale value
function result = rgb_to_grayscale(rgb_n_by_3)
    result = rgb_n_by_3 * [0.2126; 0.7152; 0.0722];
end