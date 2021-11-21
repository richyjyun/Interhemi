function [accdata, triggers, lefttrials, righttrials, fs, lefttrialsuccess, righttrialsuccess] = LoadTrainKato(fname,dwn)

%%%%%%%
% read train file

% % Reading cfg file to get fs

if (nargin<=2)
    dwn = 1;
end

train_threshold = 3000; % threshold for trigger pulse

% NEED TO CHANGE THIS TO MAKE UNIVERSAL FOR F32 FILES AS WELL
fid = fopen(fname, 'r');
if(strcmp(fname(end-2:end),'f32'))
    cfgfile = [fname(1:end-4) '.cfg']; % Change if f32 file is possibi
    load(cfgfile,'-mat');
    fs = UI.samprate;
    nchans = 9;
    skip = nchans*(dwn-1)*4;
    precision = sprintf('%d*single',nchans);
    fdata = fread(fid, [1,inf],precision,skip); %[9, inf], 'single')'; % note the transpose to put data in columns
    fdata = reshape(fdata,[nchans,length(fdata)/nchans])';
    fs = fs/dwn;
else
    fdata = fread(fid, [9, inf], 'int16=>single')'; % note the transpose to put data in columns
    fs = 1000;
end
fclose(fid);

cursorintarget = bitand(uint16(fdata(:, 9)), 1); 
targetonscreen = bitand(uint16(fdata(:, 9)), 2);
intrial = bitand(uint16(fdata(:, 9)), 4);
feeder = bitand(uint16(fdata(:, 9)), 8);
correctresponse = bitand(uint16(fdata(:, 9)), 32);
feedbackoff=bitand( uint16(fdata(:,9)), 512)/512;
boxid = floor(fdata(:, 9) / 1024); % nothing 0 reward target down 1 reward target up 2. When this goes from 2 to something else might be a good place to figure out when monkeys accomplished task.

triggers = find(fdata(:,9) > train_threshold); % andrew
intervals2 = diff([0; triggers]); % andrew
triggers = triggers(intervals2 > 2)-1; % cut out doublets 

%binarydata = cat(2, cursorintarget, targetonscreen,intrial,feeder,correctresponse,feedbackoff,boxid);

%%%%%%%%%% Calculate movement onsets
if length(unique(boxid))<4, error('not bimanual task'); end % this isnt a bimanual task (bc most likely)

targeton = boxid==2; % left trials
lefttrials = u.findOnsetsAndOffsets(targeton); % these are the trials

targeton = boxid==3; % right trials
righttrials = u.findOnsetsAndOffsets(targeton); % these are the trials

%calculate which didnt time out
tfcorrectonsets = double([0; diff(correctresponse(:)>1)>=1;]); % logical correct trial
wherecorrect = find(tfcorrectonsets)-1; % indeces where correct responses were dealt
righttrialsuccess = ismember(righttrials(:,2), wherecorrect); % right trials that are correct (where correct that dont align are experimentor dealt)
lefttrialsuccess = ismember(lefttrials(:,2), wherecorrect); % left trials that are correct

accdata = [fdata(:,1), fdata(:,4)];  %left is 1, right is 4

% % make up accel data and downsample
% downsampleby = round(fs/fs_ds);
% newfs = fs/downsampleby;
% accdata(:,1) = decimate(fdata(:,1), downsampleby);
% accdata(:,2) = decimate(fdata(:,4), downsampleby);
% 
% % get everything into ms
% 
% triggers = 1000*triggers/fs;
% lefttrials = 1000*lefttrials/fs;
% righttrials = 1000*righttrials/fs;