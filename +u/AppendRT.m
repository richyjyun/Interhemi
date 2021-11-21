function SL = AppendRT(SL, plotit)
%
% For each SL index, adds two fields called rt_l and rt_r that are as long
% as left_trials and right_trials, respectively
%
% Checking from start of trial rather than end of trial!
%
% rjy Oct 2017

if ~exist('plotit', 'var'), plotit = 0; end

window = 0.6; % sec

for i = 1:numel(SL)

    D = SL(i).Date;
    Session = char(D);
    disp(['Reaction Time Session ',Session])
    
    RFilter = u.FilterAccForMovementTimes(SL(i).accel_raw_r, SL(i).fs, 'richardson');
    LFilter = u.FilterAccForMovementTimes(SL(i).accel_raw_l, SL(i).fs, 'richardson');
   
    if isempty(RFilter)
        SL(i).rts_l = nan(size(SL(i).lefttrials,1),1);
        SL(i).rts_r = nan(size(SL(i).righttrials,1),1);
        continue;
    end
    
    LFilter = LFilter./SL(i).Max(1);
    RFilter = RFilter./SL(i).Max(2);
    threshold = 1/6; %sixth of max
    
    inds = -window*SL(i).fs:1:0.1*SL(i).fs; % indeces to look backwards and forward
    
    accelfs = 1000;
    if(isfield(SL(i),'accelfs'))
        accelfs = SL(i).accelfs;
    end
        
    % construct left successes
    leftinds = repmat(SL(i).fs/accelfs*SL(i).lefttrials(:,2)', length(inds), 1) + repmat(inds(:), 1, size(SL(i).lefttrials,1));
    % construct right successes
    rightinds = repmat(SL(i).fs/accelfs*SL(i).righttrials(:,2)', length(inds), 1) + repmat(inds(:), 1, size(SL(i).righttrials,1));
    
    % round to integers
    leftinds = floor(leftinds);
    rightinds = floor(rightinds);
    
    % remove any windows that extend before recording
    badleft = sum(leftinds<=0,1)>0;
    badright = sum(rightinds<=0,1)>0;
    
    % remove any windows that extend after recording
    badleft = badleft | leftinds(end,:)>length(LFilter);
    badright = badright | rightinds(end,:)>length(LFilter);
    
    leftinds(:,badleft) = 1;
    rightinds(:,badright) = 1;
    
    ts = linspace(-window*accelfs,0.1*accelfs, length(inds));
    
    left_move_inds = u.findfirstineachcolumn_start(LFilter(leftinds)>threshold); % for each column find first threshold crossing
    right_move_inds = u.findfirstineachcolumn_start(RFilter(rightinds)>threshold);
    
    badleft = badleft | isnan(left_move_inds); % logical vector where trial is bad
    badright = badright | isnan(right_move_inds);
    
    ts_left = nan(size(left_move_inds));
    ts_right = nan(size(right_move_inds));
    
    ts_left(~badleft) = ts(left_move_inds(~badleft)); % time relative to each trial end
    ts_right(~badright) = ts(right_move_inds(~badright));
    
    SL(i).rts_l = diff(SL(i).lefttrials, 1, 2) + ts_left(:);
    SL(i).rts_r = diff(SL(i).righttrials, 1, 2) + ts_right(:);
    
    SL(i).rts_l(SL(i).rts_l>600) = NaN; % remove those over 600ms
    SL(i).rts_r(SL(i).rts_r>600) = NaN;
    SL(i).rts_l(SL(i).rts_l<40) = NaN; % remove those under 40ms
    SL(i).rts_r(SL(i).rts_r<40) = NaN;
    
    if any(SL(i).rts_l<=0) | any(SL(i).rts_r<=0), keyboard; end % for debugging
end


if plotit
    
    packetname = 'posterfigures.ps';
    
    fs = SL(i).fs; % Hz
    
    windowsize = [-.2 10] * fs; % samples, take a wide window to see more
    
    epochs = fs*vertcat(SL(i).lefttrials, SL(i).righttrials)/1000;
    
    acceldata = horzcat(SL(i).accel_raw_l,SL(i).accel_raw_r);
    triggersignals = horzcat(LFilter(:), RFilter(:));
    
    plottingwindows = cat(2, epochs(:,1) + windowsize(1), epochs(:,1) + windowsize(2)); % windows based around target presentation
    
    %normalize to make plotting easier
    acceldata = acceldata-repmat(mean(acceldata,1),size(acceldata,1),1);
    acceldata = acceldata ./ repmat(max(acceldata,[],1), size(acceldata,1),1);
    thresh = threshold ./ max(triggersignals,[],1);
    triggersignals = triggersignals ./ repmat(max(triggersignals,[],1), size(triggersignals,1),1);
    
    epochsinc = floor(size(plottingwindows,1)/60); % sample across the whole recording
    
    figure('paperorientation', 'landscape', 'paperposition', [-1 -.7 12.5 9.6], 'position', [0 0 1100 850])
    
    subplotinds = repmat([1 3 5], 1, 40);
    
    for j = 0:39 % iterate through the plots
        
        samples = plottingwindows(5+j*epochsinc,:); % window
        
        %build the patches
        XS = cat(3, epochs,epochs);
        XS = permute(XS,[3 2 1]);
        XS = reshape(XS, 4, numel(XS)/4,1);
        YS = repmat([-1; 1; 1; -1], 1, size(XS,2));
        
        subaxis(6, 1, subplotinds(j+1), 'padding', .01, 'spacing', .01)
        patch(XS, YS, [.9 .9 .9], 'edgecolor', 'none'), hold on
        line(samples(1):samples(2), acceldata(samples(1):samples(2),1), 'linewidth', 1, 'color', [.4 .4 .4]), hold on
        line(samples(1):samples(2), triggersignals(samples(1):samples(2),1), 'linewidth', 1.5, 'color', [0 0 0]), hold on
        line([samples(1), samples(2)], [thresh(1) thresh(1)], 'color', 'r', 'linestyle', ':')
        %   line(samples(1):samples(2), double(Col9Data(samples(1):samples(2), 7))/4, 'linewidth', 1, 'color', 'b')
        
        XS = fs*SL(i).rts_l/1000;
        XS = [XS, XS]';
        YS = repmat([-1;1], 1, size(XS,2));
        
        line(XS, YS, 'linestyle', ':', 'linewidth', .7, 'color', 'r')
        
        XS = fs*SL(i).triggers/1000;
        XS = [XS, XS]';
        YS = repmat([-1;1], 1, size(XS,2));
        
        line(XS, YS, 'linestyle', ':', 'linewidth', .7, 'color', 'g')
        
        xlim(samples)
        
        
        axis off
        
        % plot the other side now
        %build the patches
        XS = cat(3, epochs,epochs);
        XS = permute(XS,[3 2 1]);
        XS = reshape(XS, 4, numel(XS)/4,1);
        YS = repmat([-1; 1; 1; -1], 1, size(XS,2));
        
        subaxis(6, 1, subplotinds(j+1)+1,'padding', .01, 'spacing', .01)
        patch(XS, YS, [.9 .9 .9], 'edgecolor', 'none'), hold on
        line(samples(1):samples(2), acceldata(samples(1):samples(2),2), 'linewidth', 1, 'color', [.4 .4 .4]), hold on
        line(samples(1):samples(2), triggersignals(samples(1):samples(2),2), 'linewidth', 1.5, 'color', [0 0 0]), hold on
        line([samples(1), samples(2)], [thresh(2) thresh(2)], 'color', 'r', 'linestyle', ':')
        
        XS = fs*SL(i).rts_r/1000;
        XS = [XS, XS]';
        YS = repmat([-1;1], 1, size(XS,2));
        
        line(XS, YS, 'linestyle', ':', 'linewidth', .7, 'color', 'r')
        
        XS = fs*SL(i).triggers/1000;
        XS = [XS, XS]';
        YS = repmat([-1;1], 1, size(XS,2));
        
        line(XS, YS, 'linestyle', ':', 'linewidth', .7, 'color', 'g')
        
        xlim(samples)
        
        axis off
        
        if mod(j+1,3)==0
            print(gcf, '-dpsc2', packetname, '-append')
            %print(gcf, '-dpdf', [packetname '.pdf'], '-append')
            clf
        end
        
    end
    
    print(gcf, '-dpsc2', packetname, '-append')
    
end