user = 'Irene';
load(['C:\Users\',user,'\Dropbox\Fetz Lab\RT\MetaData\MetaDataUbiFinal.mat'])
tempSL = SL;

for s = 1:length(tempSL)
    try
        %% Choose session and check
        SL = tempSL(s);
        fprintf('%s(%d/%d)',SL.Date,s,length(tempSL));
        
        savename = ['C:\Users\',user,'\Dropbox\Fetz Lab\RT\Ubi_Grangers_Freq\',SL.Date,'.mat'];
        
        if exist(savename)
            fprintf(' Skipped (File Exists)\n');
            continue;
        elseif isempty(SL.trig1)
            fprintf(' Skipped (No stim triggers)\n');
            continue;
        elseif SL.Bad
            fprintf(' Bad day\n');
            continue;
        end
        
        %% Load Data
        gugpath = 'R:\Yun\DualAccel\Ubi\Ubi_Guger_Data';
        trainpath = 'R:\Yun\DualAccel\Ubi\Ubi_Synced_Train';
        
        D = SL.Date;
        S = SL.Session_Guger_Train;
        Session = [char(D),'_',char(S(2))];
        trainfile = fullfile(trainpath,[Session,'.f32']);
        gugfile = fullfile(gugpath,Session);
        if(~exist(trainfile) || ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
            fprintf(' Skipped (No guger or train file)\n');
            continue;
        end
        
        fprintf(' Loading Data... '); tic;
        
        % down sampling rate
        dwn = 5;
        [accel, trig1, ~, ~, ~, ~, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi(trainfile,dwn);
        [data, fs, chnm, ~] = u.LoadGug(gugfile, dwn);
        stimstart = round(SL.trig1(1)*fs/SL.fs); stimend = round(SL.trig1(end)*fs/SL.fs);
        
        data = bpfilt(data,[1,400],fs,2);
        
        time = toc;
        fprintf('%5.2fs; ',time);
        
        %% Find channels with beta bump
        chns = 1:35; chns([17,34,35]) = []; good = [];
        [pxx,f] = pwelch(data(stimend:end,chns),[],[],[],fs);
        dip = [find(f>17,1),find(f>17,1)]; bump = [find(f>19,1),find(f>21,1)];
        for i = 1:length(chns)
            temp = smooth(pxx(:,i),1000).*f;
            temp = downsample(temp,100);
            tf = downsample(f,100);
            temp = smooth(temp,20);
            
            pks = findpeaks(temp);
            if(any(tf(pks.loc) > 15 & tf(pks.loc) < 30))
                good(end+1) = chns(i);
            end
            %         temp = smooth((smooth(pxx(:,i),1000)),1000).*f;
            %         if(mean(temp(dip(1):dip(2))) < mean(temp(bump(1):bump(2))))
            %             good(end+1) = chns(i);
            %         end
        end
        
        % use standard deviation of data to remove bad channels that snuck in
        chnstd1 = std(data(10000:(round(length(data)/10)+10000),:));
        chnstd2 = std(data((end-(round(length(data)/10)+10000)):(end-10000),:));
        good(chnstd1(good) > 100 | chnstd1(good) < 5) = [];
        good(chnstd2(good) > 100 | chnstd2(good) < 5) = [];
        
        %% Initialize some parameters
        window = [-2,2]; range = round(window(1)*fs):round(window(2)*fs);
        alpha = [8,12]; beta = [15,30]; lgamma = [65,80]; hgamma = [100,140];
        %     alpha = [30,50]; beta = [50,80]; lgamma = [30,80]; hgamma = [80,120];
        
        %% Get trial triggers (aligned with trial start)
        lefttrials = round(SL.lefttrials(:,1)*fs/SL.fs);
        lefttrials = lefttrials + round(SL.rts_l*fs/1000);
        bad = isnan(lefttrials) | ~lefttrialsuccess;
        lefttrials(bad) = [];
        
        righttrials = round(SL.righttrials*fs/SL.fs);
        righttrials(:,1) = righttrials(:,1) + round(SL.rts_r*fs/1000);
        bad = isnan(righttrials(:,1)) | ~righttrialsuccess;
        righttrials(bad,:) = [];
        
        leftpre = lefttrials(lefttrials < stimstart,1);
        leftpre = leftpre(25:end);
        leftpost = lefttrials(lefttrials > stimend,1);
        leftcond = lefttrials(lefttrials >= stimstart & lefttrials <= stimend,1);
        
        rt_l = [nanmean(SL.rts_l(lefttrials < stimstart)),...
            nanmean(SL.rts_l(lefttrials > stimend)),...
            nanmean(SL.rts_l(lefttrials >= stimstart & lefttrials <= stimend))];
        
        rightpre = righttrials(righttrials(:,1) < stimstart,1);
        rightpre = rightpre(25:end);
        rightpost = righttrials(righttrials(:,1) > stimend,1);
        rightcond = righttrials(righttrials(:,1) >= stimstart & righttrials(:,1) <= stimend,:);
        trig = SL.trig1.*fs./SL.fs;
        bad = zeros(1,length(trig));
        for i = 1:length(trig)
            temp = abs(trig(i)-rightcond(:,1));
            bad(i) = find(temp==min(temp));
        end
        rightcond(bad,:) = [];
        rightcond = rightcond(:,1);
        
        rt_r = [nanmean(SL.rts_r(righttrials(:,1) < stimstart)),...
            nanmean(SL.rts_r(righttrials(:,1) > stimend)),...
            nanmean(SL.rts_r(righttrials(:,1) >= stimstart & righttrials(:,1) <= stimend))];
        
        %% Split channels to left and right, make sure there is one in each
        lchn = good(good>16); rchn = good(good<=16);
        if(isempty(rchn))
            fprintf('Skipped (No good RH channels)\n');
            continue;
        end
        if(isempty(lchn))
            fprintf('Skipped (No good LH channels)\n');
            continue;
        end
        
        %% GC
        % X input dimensions are - channels, time, trials
        fprintf('Analyzing... ');
        tic;
        
        params.regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default) tried LWR, takes much longer but gives same results
        params.icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
        
        params.morder    = 5;%'AIC'; % 50ms, from PNAS paper
        params.momax     = 100;     % maximum model order for model order estimation
        
        params.acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)
        
        params.tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
        params.alpha_sig     = 0.005;   % significance level for significance test
        params.mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
        
        params.fs        = fs;    % sample rate (Hz)
        params.fres      = [];     % frequency resolution (empty for automatic calculation)
        
        params.nvars     = length(lchn)+length(rchn);     % number of variables. only doing pairwise
        
        params.frange    = []; %range of frequencies to look at
        
        
        gcrange = [-200,200]/1000;
        inds = find(gcrange(1)<=range/fs,1):find(gcrange(2)<=range/fs,1);
        
        params.nobs      = length(inds);   % number of observations per trial
        
        params.ntrials   = 10;
        
        X = [];
        for chn1 = [lchn,rchn]
            d = data(:,chn1);
            trialinds = getTrialinds(leftpre,range,length(data));
            trials1 = d(trialinds);
            trials1 = trials1(inds,:);
            X(end+1,:,:) = trials1;
        end
        [gc_l,pval_l] = a.GC_Final(X,params);
        
        X = [];
        for chn1 = [lchn,rchn]
            d = data(:,chn1);
            trialinds = getTrialinds(rightpre,range,length(data));
            trials1 = d(trialinds);
            trials1 = trials1(inds,:);
            X(end+1,:,:) = trials1;
        end
        [gc_r,pval_r] = a.GC_Final(X,params);
        
        time = toc;
        fprintf('%5.2fs. ',time);
        
        %% Save
        fprintf('Saving...'); tic;
        date = SL.Date;
        save(savename,'lchn','rchn','gc_l','pval_l','gc_r','pval_r','fs','-v7.3');
        
    catch
    end
    
    %% Maintenance
    time = toc;
    fprintf('%5.2fs\n',time);
    clearvars -except fname tempSL s user
    
end

