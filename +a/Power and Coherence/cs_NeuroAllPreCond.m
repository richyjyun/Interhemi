%% Save power and coherence for all frequency bands during all Experiments

user = getenv('username');
load(['C:\Users\',user,'\Dropbox\Fetz Lab\RT\MetaData\MetaDataUbiFinal.mat'])
tempSL = SL;

for s = 1:length(tempSL)
    %% Choose session and check
    SL = tempSL(s);
    fprintf('%s(%d/%d)',SL.Date,s,length(tempSL));
    
    savename = ['C:\Users\',user,'\PowerCoherence\',SL.Date,'.mat'];

    if exist(savename)
        fprintf(' Skipped (File Exists)\n');
        continue;
    elseif strcmp(SL.Notes(1:4),'100%')
        fprintf(' Skipped (100%%stim)\n');
        continue;
    elseif isempty(SL.trig1)
        fprintf(' Skipped (No stim triggers)\n');
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
    
    data = bpfilt(data,[1,500],fs,2);
    
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
    end
    
    % use standard deviation of data to remove bad channels that snuck in
    chnstd1 = std(data(10000:(round(length(data)/10)+10000),:));
    chnstd2 = std(data((end-(round(length(data)/10)+10000)):(end-10000),:));
    good(chnstd1(good) > 100 | chnstd1(good) < 5) = [];
    good(chnstd2(good) > 100 | chnstd2(good) < 5) = [];
        
    %% Initialize some parameters
    window = [-2,2]; range = round(window(1)*fs):round(window(2)*fs);
    alpha = [8,12]; beta = [15,30]; lgamma = [65,80]; hgamma = [100,140];
    
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
    
        
    %% Check and remove trials that have stimulation
    for chn = [lchn(1),rchn(1)]
        d = data(:,chn);
        trialinds = getTrialinds(leftpre,range,length(data));
        leftpredata = d(trialinds); leftpredata = leftpredata-mean(leftpredata);
        bad = any(leftpredata>500); leftpre(bad) = [];
        trialinds = getTrialinds(leftcond,range,length(data));
        leftconddata = d(trialinds); leftconddata = leftconddata-mean(leftconddata);
        bad = any(leftconddata>500); leftcond(bad) = [];
        trialinds = getTrialinds(leftpost,range,length(data));
        leftpostdata = d(trialinds); leftpostdata = leftpostdata-mean(leftpostdata);
        bad = any(leftpostdata>500); leftpost(bad) = [];
        trialinds = getTrialinds(rightpre,range,length(data));
        rightpredata = d(trialinds); rightpredata = rightpredata-mean(rightpredata);
        bad = any(rightpredata>500); rightpre(bad) = [];
        trialinds = getTrialinds(rightcond,range,length(data));
        rightconddata = d(trialinds); rightconddata = rightconddata-mean(rightconddata);
        bad = any(rightconddata>500); rightcond(bad) = [];
        trialinds = getTrialinds(rightpost,range,length(data));
        rightpostdata = d(trialinds); rightpostdata = rightpostdata-mean(rightpostdata);
        bad = any(rightpostdata>500); rightpost(bad) = [];
    end
    
    fprintf('Analyzing... ');
    tic;
    %% Power
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pow=struct;
    pow.alpha=struct; pow.beta=struct; pow.lgamma=struct; pow.hgamma=struct;
    pow.alpha.LTLH={}; pow.alpha.LTRH={}; pow.alpha.RTLH={}; pow.alpha.RTRH={};
    pow.beta.LTLH={}; pow.beta.LTRH={}; pow.beta.RTLH={}; pow.beta.RTRH={};
    pow.lgamma.LTLH={}; pow.lgamma.LTRH={}; pow.lgamma.RTLH={}; pow.lgamma.RTRH={};
    pow.hgamma.LTLH={}; pow.hgamma.LTRH={}; pow.hgamma.RTLH={}; pow.hgamma.RTRH={};
    
    erp=struct;
    erp.lefttrial=struct; erp.righttrial=struct;
    erp.lefttrial.lefthemi={}; erp.lefttrial.righthemi={}; erp.righttrial.lefthemi={}; erp.righttrial.righthemi={};
    
    fprintf('Power and ERP...');
    for chn = good
        d = data(:,chn);
        
        %% Get trialinds
        trialinds = getTrialinds(leftpre,range,length(data));
        leftpredata = d(trialinds); leftpredata = leftpredata-mean(leftpredata);
        trialinds = getTrialinds(leftcond,range,length(data));
        leftconddata = d(trialinds); leftconddata = leftconddata-mean(leftconddata);
        trialinds = getTrialinds(leftpost,range,length(data));
        leftpostdata = d(trialinds); leftpostdata = leftpostdata-mean(leftpostdata);
        trialinds = getTrialinds(rightpre,range,length(data));
        rightpredata = d(trialinds); rightpredata = rightpredata-mean(rightpredata);
        trialinds = getTrialinds(rightcond,range,length(data));
        rightconddata = d(trialinds); rightconddata = rightconddata-mean(rightconddata);
        trialinds = getTrialinds(rightpost,range,length(data));
        rightpostdata = d(trialinds); rightpostdata = rightpostdata-mean(rightpostdata);
        
        
        %% Left hemisphere
        if chn > 16
            
            sz = size(pow.alpha.LTLH,1);        
    
            % Pre/cond/post
            for e = 1:3
                if e==1
                    lefttrig = leftpredata;
                    righttrig = rightpredata;
                    erp.lefttrial.lefthemi{sz+1,e} = mean(lefttrig,2);
                    erp.righttrial.lefthemi{sz+1,e} = mean(righttrig,2);
                elseif e==2
                    lefttrig = leftconddata;
                    righttrig = rightconddata;
                    erp.lefttrial.lefthemi{sz+1,e} = mean(lefttrig,2);
                    erp.righttrial.lefthemi{sz+1,e} = mean(righttrig,2);
                else
                    lefttrig = leftpostdata;
                    righttrig = rightpostdata;
                    erp.lefttrial.lefthemi{sz+1,e} = mean(lefttrig,2);
                    erp.righttrial.lefthemi{sz+1,e} = mean(righttrig,2);
                end
                
                % Left trial
                pow.alpha.LTLH{sz+1,e} = a.getPower(lefttrig,alpha,fs);
                pow.beta.LTLH{sz+1,e} = a.getPower(lefttrig,beta,fs);
                pow.lgamma.LTLH{sz+1,e} = a.getPower(lefttrig,lgamma,fs);
                pow.hgamma.LTLH{sz+1,e} = a.getPower(lefttrig,hgamma,fs);
                
                % Right trial
                pow.alpha.RTLH{sz+1,e} = a.getPower(righttrig,alpha,fs);
                pow.beta.RTLH{sz+1,e} = a.getPower(righttrig,beta,fs);
                pow.lgamma.RTLH{sz+1,e} = a.getPower(righttrig,lgamma,fs);
                pow.hgamma.RTLH{sz+1,e} = a.getPower(righttrig,hgamma,fs);
                
            end
            
            %% Right hemisphere
        else
            
            sz = size(pow.alpha.LTRH,1);
            
            % Pre/cond/post
            for e = 1:3
                if e==1
                    lefttrig = leftpredata;
                    righttrig = rightpredata;
                    erp.lefttrial.righthemi{sz+1,e} = mean(lefttrig,2);
                    erp.righttrial.righthemi{sz+1,e} = mean(righttrig,2);
                elseif e==2
                    lefttrig = leftconddata;
                    righttrig = rightconddata;
                    erp.lefttrial.righthemi{sz+1,e} = mean(lefttrig,2);
                    erp.righttrial.righthemi{sz+1,e} = mean(righttrig,2);
                else
                    lefttrig = leftpostdata;
                    righttrig = rightpostdata;
                    erp.lefttrial.righthemi{sz+1,e} = mean(lefttrig,2);
                    erp.righttrial.righthemi{sz+1,e} = mean(righttrig,2);
                end
                
                pow.alpha.LTRH{sz+1,e} = a.getPower(lefttrig,alpha,fs);
                pow.beta.LTRH{sz+1,e} = a.getPower(lefttrig,beta,fs);
                pow.lgamma.LTRH{sz+1,e} = a.getPower(lefttrig,lgamma,fs);
                pow.hgamma.LTRH{sz+1,e} = a.getPower(lefttrig,hgamma,fs);
                
                pow.alpha.RTRH{sz+1,e} = a.getPower(righttrig,alpha,fs);
                pow.beta.RTRH{sz+1,e} = a.getPower(righttrig,beta,fs);
                pow.lgamma.RTRH{sz+1,e} = a.getPower(righttrig,lgamma,fs);
                pow.hgamma.RTRH{sz+1,e} = a.getPower(righttrig,hgamma,fs);
                
            end
            
        end
        
    end
    time = toc;
    fprintf('%5.2fs. ',time);
    
    fprintf('Coherence...');
    tic;
    %% Coherence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.tapers = [3,5]; params.Fs = fs; params.fpass = [3,120]; params.trialave = 1;
    movingwin = [0.5,0.025];
    coh=struct;
    coh.alpha=struct; coh.beta=struct; coh.lgamma=struct; coh.hgamma=struct;
    coh.alpha.LT={}; coh.alpha.RT={}; coh.beta.LT={}; coh.beta.RT={};
    coh.lgamma.LT={}; coh.lgamma.RT={}; coh.hgamma.LT={}; coh.hgamma.RT={};
    CSpect = struct; CSpect.LT = {}; CSpect.RT = {};
    for r = 1:length(rchn)
        for l = 1:length(lchn)
            
            ldata = data(:,lchn(l));
            rdata = data(:,rchn(r));
            
            %% getting trials
            trialinds = getTrialinds(leftpre,range,length(data));
            leftpreRhemi = rdata(trialinds); leftpreRhemi = leftpreRhemi-mean(leftpreRhemi);
            leftpreLhemi = ldata(trialinds); leftpreLhemi = leftpreLhemi-mean(leftpreLhemi);
            
            trialinds = getTrialinds(leftcond,range,length(data));
            leftcondRhemi = rdata(trialinds); leftcondRhemi = leftcondRhemi-mean(leftcondRhemi);
            leftcondLhemi = ldata(trialinds); leftcondLhemi = leftcondLhemi-mean(leftcondLhemi);
            
            trialinds = getTrialinds(leftpost,range,length(data));
            leftpostRhemi = rdata(trialinds); leftpostRhemi = leftpostRhemi-mean(leftpostRhemi);
            leftpostLhemi = ldata(trialinds); leftpostLhemi = leftpostLhemi-mean(leftpostLhemi);
            
            trialinds = getTrialinds(rightpre,range,length(data));
            rightpreRhemi = rdata(trialinds); rightpreRhemi = rightpreRhemi-mean(rightpreRhemi);
            rightpreLhemi = ldata(trialinds); rightpreLhemi = rightpreLhemi-mean(rightpreLhemi);
            
            trialinds = getTrialinds(rightcond,range,length(data));
            rightcondRhemi = rdata(trialinds); rightcondRhemi = rightcondRhemi-mean(rightcondRhemi);
            rightcondLhemi = ldata(trialinds); rightcondLhemi = rightcondLhemi-mean(rightcondLhemi);
            
            trialinds = getTrialinds(rightpost,range,length(data));
            rightpostRhemi = rdata(trialinds); rightpostRhemi = rightpostRhemi-mean(rightpostRhemi);
            rightpostLhemi = ldata(trialinds); rightpostLhemi = rightpostLhemi-mean(rightpostLhemi);
            
            %% Pre/Cond/Post
            
            sz = size(coh.alpha.LT,1);
            
            for e = 1:3
                if e==1
                    leftL = leftpreLhemi; leftR = leftpreRhemi;
                    rightL = rightpreLhemi; rightR = rightpreRhemi;
                elseif e==2
                    leftL = leftcondLhemi; leftR = leftcondRhemi;
                    rightL = rightcondLhemi; rightR = rightcondRhemi;
                else
                    leftL = leftpostLhemi; leftR = leftpostRhemi;
                    rightL = rightpostLhemi; rightR = rightpostRhemi;
                end
                
                % Left trial
                [C,phi,S12,S1,S2,t,f1]=cohgramc(leftL,leftR,movingwin,params);
                inds = [find(f1>=alpha(1),1),find(f1<=alpha(2),1,'last')];
                coh.alpha.LT{sz+1,e} = mean(C(:,inds(1):inds(2)),2);
                inds = [find(f1>=beta(1),1),find(f1<=beta(2),1,'last')];
                coh.beta.LT{sz+1,e} = mean(C(:,inds(1):inds(2)),2);
                inds = [find(f1>=lgamma(1),1),find(f1<=lgamma(2),1,'last')];
                coh.lgamma.LT{sz+1,e} = mean(C(:,inds(1):inds(2)),2);
                inds = [find(f1>=hgamma(1),1),find(f1<=hgamma(2),1,'last')];
                coh.hgamma.LT{sz+1,e} = mean(C(:,inds(1):inds(2)),2);
                CSpect.LT{sz+1,e} = C;

                %% Right trial
                [C,phi,S12,S1,S2,t,f1]=cohgramc(rightL,rightR,movingwin,params);
                inds = [find(f1>=alpha(1),1),find(f1<=alpha(2),1,'last')];
                coh.alpha.RT{sz+1,e} = mean(C(:,inds(1):inds(2)),2);
                inds = [find(f1>=beta(1),1),find(f1<=beta(2),1,'last')];
                coh.beta.RT{sz+1,e} = mean(C(:,inds(1):inds(2)),2);
                inds = [find(f1>=lgamma(1),1),find(f1<=lgamma(2),1,'last')];
                coh.lgamma.RT{sz+1,e} = mean(C(:,inds(1):inds(2)),2);
                inds = [find(f1>=hgamma(1),1),find(f1<=hgamma(2),1,'last')];
                coh.hgamma.RT{sz+1,e} = mean(C(:,inds(1):inds(2)),2);
                CSpect.RT{sz+1,e} = C;

            end
        end
    end
    CSpect.t = t; CSpect.f = f1;
    time = toc;
    fprintf('%5.2fs. ',time);
    
    %% Save extracted power and coherence metrics
    % Save power, coherence, reaction times, date
    fprintf('Saving...'); tic;
    date = SL.Date;
    save(savename,'pow','coh','erp','range','fs','t','rt_l','rt_r','date','f','pxx','chnstd1','chnstd2','good','CSpect','-v7.3');
    
    %% Maintenance
    time = toc;
    fprintf('%5.2fs\n',time);
    clearvars -except fname tempSL s user
    
end

diary off;
