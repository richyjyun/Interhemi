%% Check to see if frequency band amplitudes are correlated 
% If correlated, coherence argument is weaker

clear; close all;

user = getenv('username');

load(['C:\Users\',user,'\Dropbox\Fetz Lab\RT\MetaData\MetaDataUbiFinal.mat']);
Dates = extractfield(SL,'Date');

tempSL = SL;

s = 14;

SL = tempSL(s);

%% Load Data
gugpath = 'R:\Yun\DualAccel\Ubi\Ubi_Guger_Data';
trainpath = 'R:\Yun\DualAccel\Ubi\Ubi_Synced_Train';

D = SL.Date;
S = SL.Session_Guger_Train;
Session = [char(D),'_',char(S(2))];
trainfile = fullfile(trainpath,[Session,'.f32']);
gugfile = fullfile(gugpath,Session);

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
    
    [pks,locs] = findpeaks(temp);
    if(any(tf(locs) > 15 & tf(locs) < 30))
        good(end+1) = chns(i);
    end
end

% use standard deviation of data to remove bad channels that snuck in
chnstd1 = std(data(10000:(round(length(data)/10)+10000),:));
chnstd2 = std(data((end-(round(length(data)/10)+10000)):(end-10000),:));
good(chnstd1(good) > 100 | chnstd1(good) < 5) = [];
good(chnstd2(good) > 100 | chnstd2(good) < 5) = [];

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

%% Get correlations of each trial for each pairwise channel around time of trials
lchn = good(good>16); rchn = good(good<=16);
alpha = [8,12];
range = round(-1*fs):round(-0.5*fs);
trialinds = getTrialinds(sort([leftpre;rightpre]),range,length(data));

rho = {};
for r = 1:length(rchn)
    
    rdata = data(:,rchn(r));
    rdata = abs(hilbert(bpfilt(rdata,alpha,fs,3)));
    
    rdata = rdata(trialinds);
    
    for l = 1:length(lchn)
        
        disp([num2str(r),'-',num2str(l)]);
        
        ldata = data(:,lchn(l));
        ldata = abs(hilbert(bpfilt(ldata,alpha,fs,3)));
        
        ldata = ldata(trialinds);
        
        [rho{end+1},pval] = corr(rdata,ldata);
        
    end
end

avgrho = cellfun(@(x) mean(diag(x)),rho);

%% Get plv
lchn = good(good>16); rchn = good(good<=16);
alpha = [8,12];
trialinds = getTrialinds(sort([leftpre;rightpre]),range,length(data));

plv = {};
for r = 1:length(rchn)
    
    rdata = data(:,rchn(r));
    rdata = angle(hilbert(bpfilt(rdata,alpha,fs,3)));
    
    rdata = rdata(trialinds);
    
    for l = 1:length(lchn)
        
        disp([num2str(r),'-',num2str(l)]);
        
        ldata = data(:,lchn(l));
        ldata = angle(hilbert(bpfilt(ldata,alpha,fs,3)));
        
        ldata = ldata(trialinds);
        
        temp = circDiff(rdata,ldata);
        plv{end+1} = circMean(temp');
        
    end
end

















