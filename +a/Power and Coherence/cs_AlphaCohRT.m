%% Check if changes in alpha coherence is correlated with changes in reaction time
% Requested by reviewer from Frontiers

clear; close all;

user = getenv('username');

load(['C:\Users\',user,'\Dropbox\Fetz Lab\RT\MetaData\MetaDataUbiFinal.mat']);
Dates = extractfield(SL,'Date');


%% Get controls
path = ['C:\Users\',user,'\Dropbox\Fetz Lab\RT\UbiPowerCohRT_ControlFixed'];
files = dir([path,'\*.mat']);

C = {}; 
for file = 1:length(files)
    
    % Choose between type of CS
    fn = files(file).name; fn = fn(1:end-4);
    ind = cellfun(@(x) strcmp(x,fn),Dates);
    ind = find(ind);
    
    if(isempty(ind))
        continue;
    end
    d = str2num(SL(ind).Date);
    if(d == 20170329 || d==20170411 || d==20170518)
        continue;
    end
    
    load(fullfile(path,files(file).name));

    disp(fn)
    
    [keepLT,keepRT,lpair,rpair] = a.KeepCoherence(CSpect,good);
        
    CSpect.LT = CSpect.LT(keepLT,:);
    CSpect.RT = CSpect.RT(keepRT,:);
    CSpect.lpair = lpair;
    CSpect.rpair = rpair;
    
    C{end+1} = CSpect;
    
end

Cont = {}; Pairs = []; Pairs = [];
for i = 1:length(C)
    Cont = [Cont;C{i}.LT];
    Cont = [Cont;C{i}.RT];
    Pairs = [Pairs;C{i}.lpair];
    Pairs = [Pairs;C{i}.rpair];
end

chn1 = unique(Pairs(:,1));
chn2 = unique(Pairs(:,2));
Consolidate = {};
for i = 1:length(chn1)
    for j = 1:length(chn2)
        inds = Pairs(:,1)==chn1(i) & Pairs(:,2) == chn2(j);
        if(sum(inds)==0)
            continue;
        end
        temp = Cont(inds,:);
        med = {};
        for e = 1:3
            temp2 = []; 
            for t = 1:size(temp,1)
                temp2(:,:,end+1) = temp{t,e};
            end
            med{e} = sum(temp2,3)/size(temp,1);
            if(isempty(med{e}))
                keyboard;
            end
        end
        Consolidate(end+1,:) = med;
    end
end
Cont = Consolidate;


%% All data
RT = [];
COH = [];

%% For prep
C = {};
path = ['C:\Users\',user,'\Dropbox\Fetz Lab\RT\UbiPowerCohRT_FinalFixed'];
files = dir([path,'\*.mat']);
Dates = extractfield(SL,'Date');
LRT = [];
RRT = [];
for file = 1:length(files)
    
    % Choose between type of CS
    fn = files(file).name; fn = fn(1:end-4);
    ind = cellfun(@(x) strcmp(x,fn),Dates);
    ind = find(ind);
    
    if(isempty(ind))
        continue;
    elseif str2num(SL(ind).Stim_Delay) > 100 % ||  str2num(SL(ind).Stim_Delay) > 300
        continue;
    end
    
    d = str2num(SL(ind).Date);
    if(d == 20170328)
        continue;
    end    
    
    load(fullfile(path,files(file).name));

    disp(fn)
    
    [keepLT,keepRT,lpair,rpair] = a.KeepCoherence(CSpect,good);
        
    CSpect.LT = CSpect.LT(keepLT,:);
    CSpect.RT = CSpect.RT(keepRT,:);
    CSpect.lpair = lpair;
    CSpect.rpair = rpair;
    
    C{end+1} = CSpect;  
    
    LRT(end+1) = rt_l(2)-rt_l(1);
    RRT(end+1) = rt_r(2)-rt_r(1);
    
end
RT = [RT,RRT];

bnames = {'Alpha'};
bands = [8,12];
epochs = {};
siglim = 0.05;
epochs = [-0.2,0.2];
t = CSpect.t-2;
tinds = [find(t >= epochs(1),1), find(t>=epochs(2),1)];
alphachange = [];
for b = 1:length(C)
    
    inds = [find(CSpect.f>=bands(1),1),find(CSpect.f>=bands(2),1)];
    
    C1 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),Cont,'UniformOutput',false);
    C2 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),C{b}.RT,'UniformOutput',false);
    C3 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),C{b}.LT,'UniformOutput',false);

    A = cellfun(@(x,y) (y-x)./y, C1(:,1),C1(:,2),'uniformoutput',false);
    D1 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A);
    A = cellfun(@(x,y) (y-x)./y, C2(:,1),C2(:,2),'uniformoutput',false);
    D2 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A) - median(D1);
    A = cellfun(@(x,y) (y-x)./y, C3(:,1),C3(:,2),'uniformoutput',false);
    D3 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A) - median(D1);
    D1 = D1-median(D1);

    alphachange(b) = median(D3);

end

COH = [COH,alphachange];

%% For move
C = {};
path = ['C:\Users\',user,'\Dropbox\Fetz Lab\RT\UbiPowerCohRT_FinalFixed'];
files = dir([path,'\*.mat']);
Dates = extractfield(SL,'Date');
LRT = [];
RRT = [];
for file = 1:length(files)
    
    % Choose between type of CS
    fn = files(file).name; fn = fn(1:end-4);
    ind = cellfun(@(x) strcmp(x,fn),Dates);
    ind = find(ind);
    
    if(isempty(ind))
        continue;
    elseif str2num(SL(ind).Stim_Delay) < 100 ||  str2num(SL(ind).Stim_Delay) > 300
        continue;
    end
    
    d = str2num(SL(ind).Date);
    if(d == 20170328)
        continue;
    end    
    
    load(fullfile(path,files(file).name));

    disp(fn)
    
    [keepLT,keepRT,lpair,rpair] = a.KeepCoherence(CSpect,good);
        
    CSpect.LT = CSpect.LT(keepLT,:);
    CSpect.RT = CSpect.RT(keepRT,:);
    CSpect.lpair = lpair;
    CSpect.rpair = rpair;
    
    C{end+1} = CSpect;  
    
    LRT(end+1) = rt_l(2)-rt_l(1);
    RRT(end+1) = rt_r(2)-rt_r(1);
    
end
RT = [RT,LRT];

bnames = {'Alpha'};
bands = [8,12];
epochs = {};
siglim = 0.05;
epochs = [-0.2,0.2];
t = CSpect.t-2;
tinds = [find(t >= epochs(1),1), find(t>=epochs(2),1)];
alphachange = [];
for b = 1:length(C)
    
    inds = [find(CSpect.f>=bands(1),1),find(CSpect.f>=bands(2),1)];
    
    C1 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),Cont,'UniformOutput',false);
    C2 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),C{b}.RT,'UniformOutput',false);
    C3 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),C{b}.LT,'UniformOutput',false);

    A = cellfun(@(x,y) (y-x)./y, C1(:,1),C1(:,2),'uniformoutput',false);
    D1 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A);
    A = cellfun(@(x,y) (y-x)./y, C2(:,1),C2(:,2),'uniformoutput',false);
    D2 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A) - median(D1);
    A = cellfun(@(x,y) (y-x)./y, C3(:,1),C3(:,2),'uniformoutput',false);
    D3 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A) - median(D1);
    D1 = D1-median(D1);

    alphachange(b) = median(D3);

end

COH = [COH,alphachange];

%% For relax
C = {};
path = ['C:\Users\',user,'\Dropbox\Fetz Lab\RT\UbiPowerCohRT_FinalFixed'];
files = dir([path,'\*.mat']);
Dates = extractfield(SL,'Date');
LRT = [];
RRT = [];
for file = 1:length(files)
    
    % Choose between type of CS
    fn = files(file).name; fn = fn(1:end-4);
    ind = cellfun(@(x) strcmp(x,fn),Dates);
    ind = find(ind);
    
    if(isempty(ind))
        continue;
    elseif str2num(SL(ind).Stim_Delay) < 300 % ||  str2num(SL(ind).Stim_Delay) > 300
        continue;
    end
    
    d = str2num(SL(ind).Date);
    if(d == 20170328)
        continue;
    end    
    
    load(fullfile(path,files(file).name));

    disp(fn)
    
    [keepLT,keepRT,lpair,rpair] = a.KeepCoherence(CSpect,good);
        
    CSpect.LT = CSpect.LT(keepLT,:);
    CSpect.RT = CSpect.RT(keepRT,:);
    CSpect.lpair = lpair;
    CSpect.rpair = rpair;
    
    C{end+1} = CSpect;  
    
    LRT(end+1) = rt_l(2)-rt_l(1);
    RRT(end+1) = rt_r(2)-rt_r(1);
    
end
RT = [RT,LRT];

bnames = {'Alpha'};
bands = [8,12];
epochs = {};
siglim = 0.05;
epochs = [-0.2,0.2];
t = CSpect.t-2;
tinds = [find(t >= epochs(1),1), find(t>=epochs(2),1)];
alphachange = [];
for b = 1:length(C)
    
    inds = [find(CSpect.f>=bands(1),1),find(CSpect.f>=bands(2),1)];
    
    C1 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),Cont,'UniformOutput',false);
    C2 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),C{b}.RT,'UniformOutput',false);
    C3 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),C{b}.LT,'UniformOutput',false);

    A = cellfun(@(x,y) (y-x)./y, C1(:,1),C1(:,2),'uniformoutput',false);
    D1 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A);
    A = cellfun(@(x,y) (y-x)./y, C2(:,1),C2(:,2),'uniformoutput',false);
    D2 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A) - median(D1);
    A = cellfun(@(x,y) (y-x)./y, C3(:,1),C3(:,2),'uniformoutput',false);
    D3 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A) - median(D1);
    D1 = D1-median(D1);

    alphachange(b) = median(D3);

end

COH = [COH,alphachange];

%% For ipsi prep
C = {};
path = ['C:\Users\',user,'\Dropbox\Fetz Lab\RT\UbiPowerCohRT_FinalFixedIpsi'];
files = dir([path,'\*.mat']);
Dates = extractfield(SL,'Date');
LRT = [];
RRT = [];
for file = 1:length(files)
    
    % Choose between type of CS
    fn = files(file).name; fn = fn(1:end-4);
    ind = cellfun(@(x) strcmp(x,fn),Dates);
    ind = find(ind);
    
    if(isempty(ind))
        continue;
    elseif str2num(SL(ind).Stim_Delay) > 100 % ||  str2num(SL(ind).Stim_Delay) > 300
        continue;
    end
    
    d = str2num(SL(ind).Date);
    if(d == 20170328)
        continue;
    end    
    
    load(fullfile(path,files(file).name));

    disp(fn)
    
    [keepLT,keepRT,lpair,rpair] = a.KeepCoherence(CSpect,good);
        
    CSpect.LT = CSpect.LT(keepLT,:);
    CSpect.RT = CSpect.RT(keepRT,:);
    CSpect.lpair = lpair;
    CSpect.rpair = rpair;
    
    C{end+1} = CSpect;  
    
    % Calculate change in RT
    preind = SL(ind).lefttrials(:,2)<SL(ind).trig1(1);
    condind = SL(ind).lefttrials(:,2)>SL(ind).trig1(1) & SL(ind).lefttrials(:,1)<SL(ind).trig1(end);
    
    LRT(end+1) = nanmean(SL(ind).rts_l(condind)) - nanmean(SL(ind).rts_l(preind)); 
    
end

RT = [RT,LRT];

bnames = {'Alpha'};
bands = [8,12];
epochs = {};
siglim = 0.05;
epochs = [-0.2,0.2];
t = CSpect.t-2;
tinds = [find(t >= epochs(1),1), find(t>=epochs(2),1)];
alphachange = [];
for b = 1:length(C)
    
    inds = [find(CSpect.f>=bands(1),1),find(CSpect.f>=bands(2),1)];
    
    C1 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),Cont,'UniformOutput',false);
    C2 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),C{b}.RT,'UniformOutput',false);
    C3 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),C{b}.LT,'UniformOutput',false);

    A = cellfun(@(x,y) (y-x)./y, C1(:,1),C1(:,2),'uniformoutput',false);
    D1 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A);
    A = cellfun(@(x,y) (y-x)./y, C2(:,1),C2(:,2),'uniformoutput',false);
    D2 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A) - median(D1);
    A = cellfun(@(x,y) (y-x)./y, C3(:,1),C3(:,2),'uniformoutput',false);
    D3 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A) - median(D1);
    D1 = D1-median(D1);

    alphachange(b) = median(D3);

end

COH = [COH,alphachange];

%% For ipsi move
C = {};
path = ['C:\Users\',user,'\Dropbox\Fetz Lab\RT\UbiPowerCohRT_FinalFixedIpsi'];
files = dir([path,'\*.mat']);
Dates = extractfield(SL,'Date');
LRT = [];
RRT = [];
for file = 1:length(files)
    
    % Choose between type of CS
    fn = files(file).name; fn = fn(1:end-4);
    ind = cellfun(@(x) strcmp(x,fn),Dates);
    ind = find(ind);
    
    if(isempty(ind))
        continue;
    elseif str2num(SL(ind).Stim_Delay) < 100 ||  str2num(SL(ind).Stim_Delay) > 300
        continue;
    end
    
    d = str2num(SL(ind).Date);
    if(d == 20170328)
        continue;
    end    
    
    load(fullfile(path,files(file).name));

    disp(fn)
    
    [keepLT,keepRT,lpair,rpair] = a.KeepCoherence(CSpect,good);
        
    CSpect.LT = CSpect.LT(keepLT,:);
    CSpect.RT = CSpect.RT(keepRT,:);
    CSpect.lpair = lpair;
    CSpect.rpair = rpair;
    
    C{end+1} = CSpect;  

    % Calculate change in RT
    preind = SL(ind).righttrials(:,2)<SL(ind).trig1(1);
    condind = SL(ind).righttrials(:,2)>SL(ind).trig1(1) & SL(ind).righttrials(:,1)<SL(ind).trig1(end);
    
    RRT(end+1) = nanmean(SL(ind).rts_r(condind)) - nanmean(SL(ind).rts_r(preind)); 
    
    
end
RT = [RT,RRT];

bnames = {'Alpha'};
bands = [8,12];
epochs = {};
siglim = 0.05;
epochs = [-0.2,0.2];
t = CSpect.t-2;
tinds = [find(t >= epochs(1),1), find(t>=epochs(2),1)];
alphachange = [];
for b = 1:length(C)
    
    inds = [find(CSpect.f>=bands(1),1),find(CSpect.f>=bands(2),1)];
    
    C1 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),Cont,'UniformOutput',false);
    C2 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),C{b}.RT,'UniformOutput',false);
    C3 = cellfun(@(x) smooth(mean(x(:,inds(1):inds(2)),2),5),C{b}.LT,'UniformOutput',false);

    A = cellfun(@(x,y) (y-x)./y, C1(:,1),C1(:,2),'uniformoutput',false);
    D1 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A);
    A = cellfun(@(x,y) (y-x)./y, C2(:,1),C2(:,2),'uniformoutput',false);
    D2 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A) - median(D1);
    A = cellfun(@(x,y) (y-x)./y, C3(:,1),C3(:,2),'uniformoutput',false);
    D3 = cellfun(@(x) mean(x(tinds(1):tinds(2))), A) - median(D1);
    D1 = D1-median(D1);

    alphachange(b) = median(D3);

end

COH = [COH,alphachange];

%% Plot
figure; scatter(-RT,COH,'k');

xlabel('\Delta RT');
ylabel('\Delta alpha coh');

xl = xlim; yl = ylim; 
hold on; plot(xl,[0,0],'k--'); hold on; plot([0,0],yl,'k--');


keep = 1:length(RT); keep([16]) = [];
[R,P] = corrcoef(-RT(keep),COH(keep))




