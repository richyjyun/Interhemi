folder = 'C:\Users\Irene\Dropbox\Fetz Lab\RT\Ubi_Grangers';
files = dir(fullfile(folder,'*.mat')); 

GC_L = nan(length(files),32,32);
GC_R = nan(length(files),32,32);
PVals_L = nan(length(files),32,32);
PVals_R = nan(length(files),32,32);
F_L =  [];
F_R = [];
for f = 1:length(files)
    load(fullfile(folder,files(f).name));
    lchn = lchn-17; rchn = rchn+16;
    
    tempgc = nan(32,32);
    tempgc([lchn,rchn],[lchn,rchn]) = gc_l;
    GC_L(f,:,:) = tempgc;
    
    tempgc = nan(32,32);
    tempgc([lchn,rchn],[lchn,rchn]) = gc_r;
    GC_R(f,:,:) = tempgc;
    
    tempp = nan(32,32);
    tempp([lchn,rchn],[lchn,rchn]) = pval_l;
    PVals_L(f,:,:) = tempp;
    
    tempp = nan(32,32);
    tempp([lchn,rchn],[lchn,rchn]) = pval_r;
    PVals_R(f,:,:) = tempp;
    
%     tempf = nan(32,32,size(f_l,3));
%     tempf([lchn,rchn],[lchn,rchn],:) = f_l;
%     F_L(f,:,:,:) = tempf;
%     
%     tempf = nan(32,32,size(f_l,3));
%     tempf([lchn,rchn],[lchn,rchn],:) = f_l;
%     F_R(f,:,:,:) = tempf;
end

% Channel Labels
labels = cell(1,32);
for i = 1:16
    labels{i} = ['L',num2str(i)];
end
for i = 1:16
    labels{i+16} = ['R',num2str(i)];
end

% Plot GC and corresponding p-values
cmax = max(max(max(nanmean(GC_L))),max(max(nanmean(GC_R))));

figure; subplot(2,2,1);
plot_pw(squeeze(nanmean(GC_L))); title('Left Trials')
xticklabels(labels); yticklabels(labels); xtickangle(90);
colorbar; caxis([0,cmax]);

subplot(2,2,3); plot_pw(squeeze(sum(PVals_L<0.05)));
xticklabels(labels); yticklabels(labels); xtickangle(90);
colorbar; caxis([0,length(files)]);

subplot(2,2,2); plot_pw(squeeze(nanmean(GC_R))); title('Right Trials');
xticklabels(labels); yticklabels(labels); xtickangle(90);
colorbar; caxis([0,cmax]);

subplot(2,2,4); plot_pw(squeeze(sum(PVals_R<0.05)));
xticklabels(labels); yticklabels(labels); xtickangle(90);
colorbar; caxis([0,length(files)]);


%% distributions
leftgc = squeeze(nanmean(GC_L));
temp1 = leftgc(1:16,1:16); temp2 = leftgc(17:32,17:32);
samehemi = [temp1(:);temp2(:)];
temp1 = leftgc(1:16,17:32); temp2 = leftgc(17:32,1:16);
diffhemi = [temp1(:);temp2(:)];
bins = linspace(min([samehemi;diffhemi]),max([samehemi;diffhemi]),100);
figure; histogram(samehemi,bins,'normalization','probability'); 
hold on; histogram(diffhemi,bins,'normalization','probability');

rightgc = squeeze(nanmean(GC_R));
temp1 = rightgc(1:16,1:16); temp2 = rightgc(17:32,17:32);
samehemi = [temp1(:);temp2(:)];
temp1 = rightgc(1:16,17:32); temp2 = rightgc(17:32,1:16);
diffhemi = [temp1(:);temp2(:)];
bins = linspace(min([samehemi;diffhemi]),max([samehemi;diffhemi]),100);
figure; histogram(samehemi,bins,'normalization','probability'); 
hold on; histogram(diffhemi,bins,'normalization','probability');

leftsig = squeeze(sum(PVals_L<0.05));
temp1 = leftsig(1:16,1:16); temp2 = leftsig(17:32,17:32);
samehemi = [temp1(:);temp2(:)];
temp1 = leftsig(1:16,17:32); temp2 = leftsig(17:32,1:16);
diffhemi = [temp1(:);temp2(:)];
bins = 0:43;
figure; histogram(samehemi,bins,'normalization','probability'); 
hold on; histogram(diffhemi,bins,'normalization','probability');

rightsig = squeeze(sum(PVals_R<0.05));
temp1 = rightsig(1:16,1:16); temp2 = rightsig(17:32,17:32);
samehemi = [temp1(:);temp2(:)];
temp1 = rightsig(1:16,17:32); temp2 = rightsig(17:32,1:16);
diffhemi = [temp1(:);temp2(:)];
bins = 0:43;
figure; histogram(samehemi,bins,'normalization','probability'); 
hold on; histogram(diffhemi,bins,'normalization','probability');

% % Plot spectral GC (may not be worth showing)
% figure; subplot(1,2,1);
% temp = squeeze(nanmean(F_L));
% plot_spw(temp,fs,[0,100]);
% 
% subplot(1,2,2);
% temp = squeeze(nanmean(F_R));
% plot_spw(temp,fs,[0,100]);

% Plot ones with good pvalues? 
% figure; subplot(1,2,1);
% temp = squeeze(nanmean(F_L));
% plot_spw(temp,fs,[0,100]);
% 
% subplot(1,2,2);
% temp = squeeze(nanmean(F_R));
% plot_spw(temp,fs,[0,100]);




