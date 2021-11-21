%% Plot example spectral power throughout each trial 

clear; close all;

user = getenv('username');

load(['C:\Users\',user,'\Dropbox\Fetz Lab\RT\MetaData\MetaDataUbiFinal.mat']);
Dates = extractfield(SL,'Date');

tempSL = SL;

s = 12

%% Choose session and check
SL = tempSL(s);
fprintf('%s(%d/%d)',SL.Date,s,length(tempSL));

% savename = ['C:\Users\',user,'\Dropbox\Fetz Lab\RT\UbiPowerCohGO\',SL.Date,'.mat'];

% if exist(savename)
%     fprintf(' Skipped (File Exists)\n');
%     continue;
% elseif strcmp(SL.Notes(1:4),'100%')
%     fprintf(' Skipped (100%%stim)\n');
%     continue;
% elseif isempty(SL.trig1)
%     fprintf(' Skipped (No stim triggers)\n');
%     continue;
% elseif ~strcmp(SL.Condition(1:6),'Contra') || strcmp(SL.Condition,'Contra_R')
%     fprintf(' Skipped (Not a contra exp or incorrect timing)\n');
%     continue;
% end

%% Load Data
gugpath = 'R:\Yun\DualAccel\Ubi\Ubi_Guger_Data';
trainpath = 'R:\Yun\DualAccel\Ubi\Ubi_Synced_Train';

D = SL.Date;
S = SL.Session_Guger_Train;

datafile = fullfile('R:\Yun\DualAccel\Ubi',[num2str(D),'.mat']);
if(exist(datafile))
    load(datafile);
else
    
    Session = [char(D),'_',char(S(2))];
    trainfile = fullfile(trainpath,[Session,'.f32']);
    gugfile = fullfile(gugpath,Session);
    % if(~exist(trainfile) || ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
    %     fprintf(' Skipped (No guger or train file)\n');
    %     continue;
    % end
    
    fprintf(' Loading Data... '); tic;
    
    % down sampling rate
    dwn = 5;
    [accel, trig1, ~, ~, ~, ~, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi(trainfile,dwn);
    [data, fs, chnm, ~] = u.LoadGug(gugfile, dwn);
    if(~isempty(SL.trig1))
        stimstart = round(SL.trig1(1)*fs/SL.fs); stimend = round(SL.trig1(end)*fs/SL.fs);
    else
        stimend = 1;
    end
    data = bpfilt(data,[1,500],fs,2);
    
    time = toc;
    fprintf('%5.2fs; ',time);
    
    save(datafile,'accel','trig1','data','fs','chnm','stimend','-v7.3');
end

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

%% Combined spectra
% Define parameters for Chronux
params.tapers = [3,5]; params.Fs = fs; params.fpass = [3,55]; params.trialave = 1;
movingwin = [0.5,0.01];

% w/respect to GO signal
trials = sort([SL.lefttrials(:,1)]);
trials = trials(trials < SL.trig1(1) | trials > SL.trig1(end));

trig = round(trials*fs/SL.fs);
range = round(-2*fs):round(2*fs);
trialinds = getTrialinds(trig,range,length(data));

Sc = []; Si = []; rh = good;
for i = 1:length(rh)
    disp(rh(i));
    filt = bpfilt(data(:,rh(i)),[.5,60],fs,2);
    trials = filt(trialinds); trials = trials-mean(trials);
    % Moving window spectrum
    if(rh(i) > 16)
        [Sc(:,:,end+1),t,f] = mtspecgramc(trials,movingwin,params);
    else
        [Si(:,:,end+1),t,f] = mtspecgramc(trials,movingwin,params);
    end
end

Scontra = mean(Sc,3);
Sipsi = mean(Si,3);

% Interpolated to smooth it out
interpS = interp2(pow2db(Scontra),5);

tempt = linspace(t(1),t(end),size(interpS,1));
tempf = linspace(f(1),f(end),size(interpS,2))-2;
figure; subplot(2,2,1);
imagesc(tempt-2,tempf,interpS');
xlim([-1,1.5]); ylim([5,50]);
set(gca,'YDir','normal'); colormap('turbo'); 
xlabel('Time from GO (s)'); ylabel('Frequency (Hz)');
c = colorbar; c.Label.Position(1) = 2; c.Label.String = 'Power (dB)';
c1 = caxis; title('Contralateral Hemisphere');

yl = ylim; hold on; plot([0,0],yl,'k--');

subplot(2,2,3);
imagesc(tempt-2,tempf,zscore(interpS)')
xlim([-1,1.5]); ylim([5,50]);
set(gca,'YDir','normal'); colormap('turbo'); 
xlabel('Time from GO (s)'); ylabel('Frequency (Hz)');
c = colorbar; c.Label.Position(1) = 2; c.Label.String = 'Normalized Power (a.u.)';
yl = ylim; hold on; plot([0,0],yl,'k--');
c2 = caxis; title('Contralateral Hemisphere');

% Interpolated to smooth it out
interpS = interp2(pow2db(Sipsi),5);

subplot(2,2,2);
imagesc(tempt-2,tempf,interpS');
xlim([-1,1.5]); ylim([5,50]);
set(gca,'YDir','normal'); colormap('turbo'); 
xlabel('Time from GO (s)'); ylabel('Frequency (Hz)');
c = colorbar; c.Label.Position(1) = 2; c.Label.String = 'Power (dB)';
caxis(c1); title('Ipsilateral Hemisphere');

yl = ylim; hold on; plot([0,0],yl,'k--');

subplot(2,2,4);
imagesc(tempt-2,tempf,zscore(interpS)')
xlim([-1,1.5]); ylim([5,50]);
set(gca,'YDir','normal'); colormap('turbo'); 
xlabel('Time from GO (s)'); ylabel('Frequency (Hz)');
c = colorbar; c.Label.Position(1) = 2; c.Label.String = 'Normalized Power (a.u.)';
yl = ylim; hold on; plot([0,0],yl,'k--');
caxis(c2); title('Ipsilateral Hemisphere');

%% Coherence
% Define parameters for Chronux
params.tapers = [3,5]; params.Fs = fs; params.fpass = [3,55]; params.trialave = 1;
movingwin = [0.5,0.01];

% w/respect to GO signal
trials = sort([SL.righttrials(:,1);SL.lefttrials(:,1)]);
trials = trials(trials < SL.trig1(1) | trials > SL.trig1(end));

trig = round(trials*fs/SL.fs);
range = round(-2*fs):round(2*fs);
trialinds = getTrialinds(trig,range,length(data));

lchns = good(good <= 16);
rchns = good(good > 16);

C = [];
for i = 13:length(lchns)
    for j = 1:length(rchns)
        disp(['Chns: ',num2str(lchns(i)),'-',num2str(rchns(j))]);
        temp = data(:,lchns(i));
        trials1 = temp(trialinds); trials1 = trials1-mean(trials1);
        temp = data(:,rchns(j));
        trials2 = temp(trialinds); trials2 = trials2-mean(trials2);
        [C(:,:,end+1),phi,S12,S1,S2,t,f] = cohgramc(trials1,trials2,movingwin,params);
    end
end
        
combC = mean(C,3);

% plot
interpC = interp2(combC,5);

tempt = linspace(t(1),t(end),size(interpC,1));
tempf = linspace(f(1),f(end),size(interpC,2))-2;
figure; subplot(1,2,1); imagesc(tempt-2,tempf,interpC')
xlim([-1,1.5]); ylim([5,50]); caxis([0,0.3]);
yl = ylim; hold on; plot([0,0],yl,'k--');
set(gca,'YDir','normal'); colormap('turbo'); 
xlabel('Time from GO (s)'); ylabel('Frequency (Hz)');
c = colorbar; c.Label.Position(1) = 2; c.Label.String = 'Coherence';

subplot(1,2,2); imagesc(tempt-2,tempf,zscore(interpC)')
xlim([-1,1.5]); ylim([5,50]);
yl = ylim; hold on; plot([0,0],yl,'k--');
set(gca,'YDir','normal'); colormap('turbo'); 
xlabel('Time from GO (s)'); ylabel('Frequency (Hz)');
c = colorbar; c.Label.Position(1) = 2; c.Label.String = 'Normalized Coherence';

datafile = ['C:\Users\Richy Yun\Dropbox\Fetz Lab\RT\Manuscript\figures\RichyFigures\Spectra\Coherence',...
    num2str(D),'.mat'];
save(datafile,'interpC','tempt','tempf','C','lchns','rchns');







