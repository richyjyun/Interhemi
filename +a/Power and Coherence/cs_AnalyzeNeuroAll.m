%% Plot aggregates of all experiments

user = getenv('username');

load(['C:\Users\',user,'\Dropbox\Fetz Lab\RT\MetaData\MetaDataUbiConditions.mat']);

path = ['C:\Users\',user,'\Dropbox\Fetz Lab\RT\UbiPowerCohRT_FinalFixed2'];
files = dir([path,'\*.mat']);

POW = {};
COH = {};
ERP = {};
rtl = [];
rtr = [];

Dates = extractfield(SL,'Date');

for f = 1:length(files)
    
    % Choose between type of CS
    fn = files(f).name; fn = fn(1:end-4);
    ind = cellfun(@(x) strcmp(x,fn),Dates);
    ind = find(ind);
    
    if(strcmp(fn, '20170405'))
        continue;
    end
    
    d = str2num(SL(ind).Date);
    if(d == 20170329 || d==20170411 || d==20170518)
        continue;
    end
    
    if(isempty(ind))
        continue;
    elseif str2num(SL(ind).Stim_Delay) < 300 %||  str2num(SL(ind).Stim_Delay) > 500
        continue;
    end
    
    load(fullfile(path,files(f).name));
    
    if(size(pow.alpha.LTLH,2) < 3 || size(coh.alpha.LT,2) < 3)
        continue;
    end
    
    disp(fn)
    
    POW{end+1} = pow;
    COH{end+1} = coh;
    ERP{end+1} = erp;
    rtl(end+1,:) = rt_l;
    rtr(end+1,:) = rt_r;
    
end

rtl = mean(rtl); rtr = mean(rtr);
xl = [-1.75,1.75];
t = t-2;

fname = ['C:\Users\',user,'\Dropbox\Fetz Lab\RT\PowerCohPackets_Fixed\Contra_Pre.ps'];

%% ERP
plotERPFig(range/fs,ERP,'lefttrial',fname);
plotERPFig(range/fs,ERP,'righttrial',fname);

%% Power - Alpha
plotPowFig(range/fs,POW,'Power','alpha',rtl,rtr,xl,fname);

%% Power - Beta
plotPowFig(range/fs,POW,'Power','beta',rtl,rtr,xl,fname);

%% Power - Low Gamma
plotPowFig(range/fs,POW,'Power','lgamma',rtl,rtr,xl,fname);

%% Power - High Gamma
plotPowFig(range/fs,POW,'Power','hgamma',rtl,rtr,xl,fname);

%% Coherence - Alpha
plotCohFig(t,COH,'Coherence','alpha',rtl,rtr,xl,fname);

%% Coherence - Beta
plotCohFig(t,COH,'Coherence','beta',rtl,rtr,xl,fname);

%% Coherence - Low Gamma
plotCohFig(t,COH,'Coherence','lgamma',rtl,rtr,xl,fname);

%% Coherence - High Gamma
plotCohFig(t,COH,'Coherence','hgamma',rtl,rtr,xl,fname);


%% Convert to pdf
callps2pdf(fname);

%% General Plotting
% Plot power
function plotPowFig(x,data,datatype,band,rt_l,rt_r,xl,fname)

figure('visible','off');
subplot(2,1,1);
singlePlot(x,data,datatype,band,'LTLH',rt_l,xl);
subplot(2,1,2); 
singlePlot(x,data,datatype,band,'LTRH',rt_l,xl);
print('-fillpage','-painters',gcf, '-dpsc2', fname, '-append');
close(gcf);

figure('visible','off');
subplot(2,1,1);
singlePlot(x,data,datatype,band,'RTLH',rt_r,xl);
subplot(2,1,2);
singlePlot(x,data,datatype,band,'RTRH',rt_r,xl);
print('-fillpage','-painters',gcf, '-dpsc2', fname, '-append');
close(gcf);

end

% Plot coherence
function plotCohFig(x,data,datatype,band,rt_l,rt_r,xl,fname)

figure('visible','off');
subplot(2,1,1);
singlePlot(x,data,datatype,band,'LT',rt_l,xl);
subplot(2,1,2);
singlePlot(x,data,datatype,band,'RT',rt_r,xl);
print('-fillpage','-painters',gcf, '-dpsc2', fname, '-append');
close(gcf);

end

% Plot ERP
function plotERPFig(x,data,trial,fname)

figure('visible','off');
subplot(2,1,1); 
temp = cellfun(@(x) x.(trial).lefthemi,data,'uniformoutput',0);
for i = 1:length(temp)
    temp2 = temp{i};
    for j = 1:3
        temp3{i,j} = mean(zscore(cell2mat(temp2(:,j)')),2);
    end
end
for i=1:size(temp3,2)
    temp4{i} = mean(cell2mat(temp3(:,i)'),2);
end
plot(x,temp4{1}); hold on; plot(x,temp4{2}); plot(x,temp4{3});
title([trial,' LeftHemi'])
xlabel('Time (s)')

subplot(2,1,2); 
temp = cellfun(@(x) x.(trial).righthemi,data,'uniformoutput',0);
for i = 1:length(temp)
    temp2 = temp{i};
    for j = 1:3
        temp3{i,j} = mean(zscore(cell2mat(temp2(:,j)')),2);
    end
end
for i=1:size(temp3,2)
    temp4{i} = mean(cell2mat(temp3(:,i)'),2);
end
plot(x,temp4{1}); hold on; plot(x,temp4{2}); plot(x,temp4{3});
title([trial,' RightHemi'])
xlabel('Time (s)')

print('-fillpage','-painters',gcf, '-dpsc2', fname, '-append');
close(gcf);

end

% Make a plot
function singlePlot(x,data,datatype,band,type,rt,xl)

avg = calcAverage(data,band,type);
if(length(avg) > 1000)
    smt = 300;
else
    smt = 5;
end
for i=1:size(avg,2)
    avg(:,i) = smoothdata(avg(:,i), 'gaussian', smt);
end
plot(x,avg); xlim(xl); yl = ylim; hold on;

xlabel('Time since movement onset (s)'); ylabel([band,' ',datatype]);
title(type);

end

% Average the values 
function avg = calcAverage(data,band,type)

% Normalize off of pre 
mu = cellfun(@(x) mean(mean(cat(2,x.(band).(type){:,1}),2)),data,'uniformoutput',0);
sigma = cellfun(@(x) std(mean(cat(2,x.(band).(type){:,1}),2)),data,'uniformoutput',0);

for i = 1:size(data{1}.(band).(type),2)
    temp = cellfun(@(x,y,z) median(cat(2,x.(band).(type){:,i}),2),data,mu,sigma,'uniformoutput',0);
    temp = (cell2mat(temp));
    avg(:,i) = smooth(mean(temp,2));
end
    
end
