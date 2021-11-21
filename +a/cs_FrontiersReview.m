%% Figures for Frontiers Reviewers

%% Distribution of stimulation relative to RT 
user = 'Richy Yun';

if(~exist('ubiSL'))
    load(['C:\Users\',user,'\Dropbox\Fetz Lab\RT\MetaData\MetaDataUbiFinal.mat']);
    ubiSL = SL;
end

if(~exist('igorSL'))
    load(['C:\Users\',user,'\Dropbox\Fetz Lab\RT\MetaData\MetaDataIgorFinal.mat']);
    igorSL = SL;
end

if(~exist('katoSL'))
    load(['C:\Users\',user,'\Dropbox\Fetz Lab\RT\MetaData\MetaDataKatoFinal.mat']);
    katoSL = SL;
end

dt = [];
for i = 1:3
    switch i
        case 1
            SL = ubiSL;
        case 2
            SL = igorSL;
        case 3
            SL = katoSL;
    end
    for s = 1:length(SL)
        if(isempty(SL(s).trig1) || strcmp(SL(s).Condition,'Control') || strcmp(SL(s).Condition,'NaN') || ~isempty(SL(s).Bad) || strcmp(SL(s).Condition(end),'R')...
                || strcmp(SL(s).Condition,'nostim'))
            continue;
        end
        
        stim = 0;
        if strcmp(SL(s).Condition(end),'M')
            bins = sort([SL(s).lefttrials(:,1);SL(s).righttrials(:,1)]);
            inds = discretize(SL(s).trig1,bins);
            stim = median(SL(s).trig1-bins(inds));
        else
            stim = str2num(SL(s).Stim_Delay);
        end
        
        stimstart = SL(s).trig1(1);
        
        left = find(SL(s).lefttrials(:,2) < stimstart,1,'last');
        right = find(SL(s).righttrials(:,2) < stimstart,1,'last');
        
        leftrt = SL(s).rts_l(25:left);
        rightrt = SL(s).rts_r(25:right);
        
        rts = [leftrt;rightrt];
        rt = nanmedian(rts);
        
        dt(end+1) = stim - rt;
                
    end
end

figure; 
scatter(dt,zeros(1,length(dt)),'k'); ylim([-0.5,1]);
yl = ylim; hold on; plot([0,0],yl,'k--'); plot([300,300],yl,'k--');
xlim([-300,600]); xlabel('Time from RT (ms)');
yticks([]);
text(-150,0.5,'CS_{prep}','horizontalalignment','center','fontsize',12);
text(150,0.5,'CS_{move}','horizontalalignment','center','fontsize',12);
text(450,0.5,'CS_{relax}','horizontalalignment','center','fontsize',12);
set(gca,'FontSize',10)
temp = gca;
temp.YAxis.Visible = 'off';

%% Checking how many trials are removed due to whiskers / outliers
RRTs = {};
LRTs = {};
for i = 1:3
    switch i
        case 1
            SL = ubiSL;
        case 2
            SL = igorSL;
        case 3
            SL = katoSL;
    end
    for s = 1:length(SL)
        if(isempty(SL(s).trig1) || strcmp(SL(s).Condition,'Control') || strcmp(SL(s).Condition,'NaN'))
            continue;
        end
        
        stimstart = SL(s).trig1(1);
        
        left = find(SL(s).lefttrials(:,2) < stimstart,1,'last');
        right = find(SL(s).righttrials(:,2) < stimstart,1,'last');
        
        leftrt = SL(s).rts_l(25:left);
        rightrt = SL(s).rts_r(25:right);
        
        LRTs{end+1} = leftrt;
        RRTs{end+1} = rightrt;
                
    end
end

RTs = RRTs;

q1 = cellfun(@(x) prctile(x,25), RTs,'uniformoutput',false);
q3 = cellfun(@(x) prctile(x,75), RTs,'uniformoutput',false);

w = 1.5;
upper = cellfun(@(x,y) y+w.*(y-x),q1,q3,'uniformoutput',false);
lower = cellfun(@(x,y) x-w.*(y-x),q1,q3,'uniformoutput',false);

prc = cellfun(@(x,y,z) sum(x<y | x>z)./length(x),RTs,lower,upper);

% Looks like roughly 5%

%% ERP analysis
%% Load SL
user = 'Richy Yun';

load(['C:\Users\',user,'\Dropbox\Fetz Lab\RT\MetaData\MetaDataUbiConditions.mat']);

path = ['C:\Users\',user,'\Dropbox\Fetz Lab\RT\UbiERPAccel'];
files = dir([path,'\*.mat']);

Dates = extractfield(SL,'Date');

%% Load Data

prepERPs = {};
moveERPs = {};
postERPs = {};

for f = 1:length(files)
    
    % Choose between type of CS
    fn = files(f).name; fn = fn(1:end-4);
    ind = cellfun(@(x) strcmp(x,fn),Dates);
    ind = find(ind);
    
    if(isempty(ind)) %ind==16 has bad LFPs
        continue;
    end
    
    load(fullfile(path,files(f).name));
    
    if str2num(SL(ind).Stim_Delay) < 100 %||  str2num(SL(ind).Stim_Delay) > 300
        prepERPs{1,end+1} = L_ERP;
        prepERPs{2,end} = R_ERP;
    elseif str2num(SL(ind).Stim_Delay) > 300
        postERPs{1,end+1} = L_ERP;
        postERPs{2,end} = R_ERP;
    else
        moveERPs{1,end+1} = L_ERP;
        moveERPs{2,end} = R_ERP;
    end
    
end

rh = 1:16;
lh = 16:32;

postERPs(:,4) = []; % bad LFPs

%% Prep
Ltrial_ERPs = {};
Rtrial_ERPs = {};
for i = 1:length(prepERPs)
    
    %% left trials
    temp = prepERPs{1,i};
    
    bad = any(squeeze(any(temp>500))') | any(squeeze(any(temp<-500))');
    bad = bad | any(all(temp<5,3),1);
    bad(33:35) = 0;
    bad = find(bad);
    
    [~,ind,~] = intersect(lh,bad);
    templh = lh;
    templh(ind) = [];
    
    % left hemisphere
    temp_l = temp(:,templh,:);
    
    left = [];
    for e = 1:3
        left(e,:) = mean(squeeze(temp_l(e,:,:)));
    end
    Ltrial_ERPs{1,i} = left;

    % right hemisphere
    [~,ind,~] = intersect(rh,bad);
    temprh = rh;
    temprh(ind) = [];
    temp_r = temp(:,temprh,:);
    
    right = [];
    for e = 1:3
        right(e,:) = mean(squeeze(temp_r(e,:,:)));
    end
    Ltrial_ERPs{2,i} = right;
    
    
    %% right trials
    temp = prepERPs{2,i};
    
    bad = any(squeeze(any(temp>500))') | any(squeeze(any(temp<-500))');
    bad = bad | any(all(temp<5,3),1);
    bad(33:35) = 0;
    bad = find(bad);
    
    [~,ind,~] = intersect(lh,bad);
    templh = lh;
    templh(ind) = [];
    
    % left hemisphere
    temp_l = temp(:,templh,:);
    
    left = [];
    for e = 1:3
        left(e,:) = mean(squeeze(temp_l(e,:,:)));
    end
    Rtrial_ERPs{1,i} = left;

    % right hemisphere
    [~,ind,~] = intersect(rh,bad);
    temprh = rh;
    temprh(ind) = [];
    temp_r = temp(:,temprh,:);
    
    right = [];
    for e = 1:3
        right(e,:) = mean(squeeze(temp_r(e,:,:)));
    end
    Rtrial_ERPs{2,i} = right;
end

% Combine
figure;
for e = 1:3
    % Left trials (ipsi trials)
    temp = cellfun(@(x) x(e,:), Ltrial_ERPs, 'UniformOutput', false);
    % Left hemisphere (ipsi hemi)
    temp_l = cell2mat(temp(1,:)');
    subplot(2,2,1); hold on; plot(range/fs,mean(temp_l));
    title('Ipsi Trials'); ylabel('Ipsi Hemisphere');
    % Right hemisphere (contra hemi)
    temp_r = cell2mat(temp(2,:)');
    subplot(2,2,3); hold on; plot(range/fs,mean(temp_r));
    ylabel('Contra Hemisphere');
    
    % Right trials (contra trials)
    temp = cellfun(@(x) x(e,:), Rtrial_ERPs, 'UniformOutput', false);
    % Left hemisphere (ipsi hemi)
    temp_l = cell2mat(temp(1,:)');
    subplot(2,2,2); hold on; plot(range/fs,mean(temp_l));
    title('Contra Trials');
    % Right hemisphere (contra hemi)
    temp_r = cell2mat(temp(2,:)');
    subplot(2,2,4); hold on; plot(range/fs,mean(temp_r));
end

%% Move
Ltrial_ERPs = {};
Rtrial_ERPs = {};
for i = 1:length(moveERPs)
    
    %% left trials
    temp = moveERPs{1,i};
    
    bad = any(squeeze(any(temp>500))') | any(squeeze(any(temp<-500))');
    bad = bad | any(all(temp<5,3),1);
    bad(33:35) = 0;
    bad = find(bad);
    
    [~,ind,~] = intersect(lh,bad);
    templh = lh;
    templh(ind) = [];
    
    % left hemisphere
    temp_l = temp(:,templh,:);
    
    left = [];
    for e = 1:3
        left(e,:) = mean(squeeze(temp_l(e,:,:)));
    end
    Ltrial_ERPs{1,i} = left;

    % right hemisphere
    [~,ind,~] = intersect(rh,bad);
    temprh = rh;
    temprh(ind) = [];
    temp_r = temp(:,temprh,:);
    
    right = [];
    for e = 1:3
        right(e,:) = mean(squeeze(temp_r(e,:,:)));
    end
    Ltrial_ERPs{2,i} = right;
    
    
    %% right trials
    temp = moveERPs{2,i};
    
    bad = any(squeeze(any(temp>500))') | any(squeeze(any(temp<-500))');
    bad = bad | any(all(temp<5,3),1);
    bad(33:35) = 0;
    bad = find(bad);
    
    [~,ind,~] = intersect(lh,bad);
    templh = lh;
    templh(ind) = [];
    
    % left hemisphere
    temp_l = temp(:,templh,:);
    
    left = [];
    for e = 1:3
        left(e,:) = mean(squeeze(temp_l(e,:,:)));
    end
    Rtrial_ERPs{1,i} = left;

    % right hemisphere
    [~,ind,~] = intersect(rh,bad);
    temprh = rh;
    temprh(ind) = [];
    temp_r = temp(:,temprh,:);
    
    right = [];
    for e = 1:3
        right(e,:) = mean(squeeze(temp_r(e,:,:)));
    end
    Rtrial_ERPs{2,i} = right;
end

% Combine
figure;
for e = 1:3
    % Left trials (ipsi trials)
    temp = cellfun(@(x) x(e,:), Ltrial_ERPs, 'UniformOutput', false);
    % Left hemisphere (ipsi hemi)
    temp_l = cell2mat(temp(1,:)');
    subplot(2,2,1); hold on; plot(range/fs,mean(temp_l));
    title('Ipsi Trials'); ylabel('Ipsi Hemisphere');
    % Right hemisphere (contra hemi)
    temp_r = cell2mat(temp(2,:)');
    subplot(2,2,3); hold on; plot(range/fs,mean(temp_r));
    ylabel('Contra Hemisphere');
    
    % Right trials (contra trials)
    temp = cellfun(@(x) x(e,:), Rtrial_ERPs, 'UniformOutput', false);
    % Left hemisphere (ipsi hemi)
    temp_l = cell2mat(temp(1,:)');
    subplot(2,2,2); hold on; plot(range/fs,mean(temp_l));
    title('Contra Trials');
    % Right hemisphere (contra hemi)
    temp_r = cell2mat(temp(2,:)');
    subplot(2,2,4); hold on; plot(range/fs,mean(temp_r));
end


%% Post
Ltrial_ERPs = {};
Rtrial_ERPs = {};
for i = 1:length(postERPs)
    
    %% left trials
    temp = postERPs{1,i};
    
    bad = any(squeeze(any(temp>500))') | any(squeeze(any(temp<-500))');
    bad = bad | any(all(temp<5,3),1);
    bad(33:35) = 0;
    bad = find(bad);
    
    [~,ind,~] = intersect(lh,bad);
    templh = lh;
    templh(ind) = [];
    
    % left hemisphere
    temp_l = temp(:,templh,:);
    
    left = [];
    for e = 1:3
        left(e,:) = mean(squeeze(temp_l(e,:,:)));
    end
    Ltrial_ERPs{1,i} = left;

    % right hemisphere
    [~,ind,~] = intersect(rh,bad);
    temprh = rh;
    temprh(ind) = [];
    temp_r = temp(:,temprh,:);
    
    right = [];
    for e = 1:3
        right(e,:) = mean(squeeze(temp_r(e,:,:)));
    end
    Ltrial_ERPs{2,i} = right;
    
    
    %% right trials
    temp = postERPs{2,i};
    
    bad = any(squeeze(any(temp>500))') | any(squeeze(any(temp<-500))');
    bad = bad | any(all(temp<5,3),1);
    bad(33:35) = 0;
    bad = find(bad);
    
    [~,ind,~] = intersect(lh,bad);
    templh = lh;
    templh(ind) = [];
    
    % left hemisphere
    temp_l = temp(:,templh,:);
    
    left = [];
    for e = 1:3
        left(e,:) = mean(squeeze(temp_l(e,:,:)));
    end
    Rtrial_ERPs{1,i} = left;

    % right hemisphere
    [~,ind,~] = intersect(rh,bad);
    temprh = rh;
    temprh(ind) = [];
    temp_r = temp(:,temprh,:);
    
    right = [];
    for e = 1:3
        right(e,:) = mean(squeeze(temp_r(e,:,:)));
    end
    Rtrial_ERPs{2,i} = right;
end

% Combine
figure;
for e = 1:3
    % Left trials (ipsi trials)
    temp = cellfun(@(x) x(e,:), Ltrial_ERPs, 'UniformOutput', false);
    % Left hemisphere (ipsi hemi)
    temp_l = cell2mat(temp(1,:)');
    subplot(2,2,1); hold on; plot(range/fs,mean(temp_l));
    title('Ipsi Trials'); ylabel('Ipsi Hemisphere');
    % Right hemisphere (contra hemi)
    temp_r = cell2mat(temp(2,:)');
    subplot(2,2,3); hold on; plot(range/fs,mean(temp_r));
    ylabel('Contra Hemisphere');
    
    % Right trials (contra trials)
    temp = cellfun(@(x) x(e,:), Rtrial_ERPs, 'UniformOutput', false);
    % Left hemisphere (ipsi hemi)
    temp_l = cell2mat(temp(1,:)');
    subplot(2,2,2); hold on; plot(range/fs,mean(temp_l));
    title('Contra Trials');
    % Right hemisphere (contra hemi)
    temp_r = cell2mat(temp(2,:)');
    subplot(2,2,4); hold on; plot(range/fs,mean(temp_r));
end








