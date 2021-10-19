function csl_out = calculate_csl(twin_map, active_variants, trace_dir)

for iTwin = 1:6
    if active_variants(iTwin)>0
        map_r = imrotate(twin_map, trace_dir(iTwin), 'nearest', 'loose');   % rotated map
        [nr,nc] = size(map_r);
        label_map = zeros(nr,nc);    % to store assigned temporary label
        csl_map = zeros(nr,nc);
        for ir = 1:nr
            if any(map_r(ir,:))
                icL_back = find(map_r(ir,:),1,'first');
                icR_back = find(map_r(ir,:),1,'last');
                icL_front = find(map_r(ir,:),1,'first');
                icR_front = find(map_r(ir,:),1,'last');
                while (icL_front<=icR_front)
                    csl_length = 0;
                    if (icL_front-icL_back)<=(icR_back-icR_front)
                        % search for connected segments from left to right
                        while(map_r(ir,icL_front))
                            label_map(ir,icL_front) = 1;  % prepare to assign label
                            map_r(ir,icL_front) = 0;     % make element on map_r 0 (temp change)
                            if icL_front < nc
                                icL_front = icL_front + 1;   % move pointer forward to the right
                                csl_length = csl_length + 1;
                            end
                        end
                        csl_map(label_map==1) = csl_length;
                        label_map(label_map==1) = 0;     % change back (equivalent to change to gb label)
                        icL_back = icL_front - 1;    % assign left side back
                        icL_front = find(map_r(ir,:),1,'first'); % search left side front again
                    else
                        while (map_r(ir,icR_front))
                            label_map(ir,icR_front) = 1;
                            map_r(ir,icR_front) = 0;
                            if icR_front > 1
                                icR_front = icR_front - 1;
                                csl_length = csl_length + 1;
                            end
                        end
                        csl_map(label_map==1) = csl_length;
                        label_map(label_map==1) = 0;
                        icR_back = icR_front + 1;
                        icR_front = find(map_r(ir,:),1,'last');
                    end
                end
            end
        end
        
        temp = imrotate(csl_map, -trace_dir(iTwin), 'nearest', 'loose');
        img1_template = (twin_map>0);
        img2_signal = (temp>0);
        [yOffSet, xOffSet] = normxcorr2A_register(img1_template, img2_signal, [0 0 0 0], [0 0 0 0], 0);
        indC_back_min = 1 + xOffSet;
        indC_back_max = indC_back_min + size(twin_map,2)-1;
        indR_back_min = 1 + yOffSet;
        indR_back_max = indR_back_min + size(twin_map,1)-1;
        
        try
            csl_p(:,:,iTwin) = temp(indR_back_min:indR_back_max, indC_back_min:indC_back_max);
        catch
            csl_p(:,:,iTwin) = zeros(size(twin_map));   % in case temp is all zero, and cannot find correlation.  
        end
    end
end

try
    [nr,nc,~] = size(csl_p);
    csl = zeros(nr,nc);   % variant number map of the pixels need to analyze
    for ir=1:nr
        for ic=1:nc
            [max_value, vN] = max(csl_p(ir,ic,:));
            if max_value>0
                csl(ir,ic) = vN;
            end
        end
    end
catch
    % if no active variants, just all 0
    csl = zeros(size(twin_map));
end

% remove those with twin_map==0
csl(twin_map==0) = 0;
% myplot(csl);
csl_out = one_pass_fill(csl);
% myplot(csl_out);
if sum(csl_out(twin_map(:)==0))
   error('need check, did not clean enough') 
end
end

