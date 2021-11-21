function [accdata, trig1, trig2, lefttrials, righttrials, fs, lefttrialsuccess, righttrialsuccess] = LoadTrainUbi(fname,dwn)
% this function loads ubi train data. Modified from LoadI16Bimanual, can
% load either i16 or f32.
%
% RJY 2/21/17

if (nargin<2)
    dwn = 1;
end

train_threshold = 3000; % threshold for trigger pulse

fid = fopen(fname, 'r');
if(strcmp(fname(end-2:end),'f32'))
    nchans = 9;
    skip = nchans*(dwn-1)*4;
    precision = sprintf('%d*single',nchans);
    fdata = fread(fid, [1,inf],precision,skip); %[9, inf], 'single')'; % note the transpose to put data in columns
    fdata = reshape(fdata,[nchans,length(fdata)/nchans])';
    fs = 9600/dwn;
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

trig1 = find(bitand(uint16(fdata(:,9)),64) > 0); % andrew
intervals2 = diff([0; trig1]); % andrew
trig1 = trig1(intervals2 > 2)-1; % cut out doublets
if(isempty(trig1))
    trig1 = find(fdata(:,7) > train_threshold); % andrew
    intervals2 = diff([0; trig1]); % andrew
    trig1 = trig1(intervals2 > 2)-1; % cut out doublets
end

trig2 = find(bitand(uint16(fdata(:,9)),128) > 0);
intervals2 = diff([0; trig2]); % andrew
trig2 = trig2(intervals2 > 2)-1; % cut out doublets

binarydata = cat(2, cursorintarget, targetonscreen,intrial,feeder,correctresponse,feedbackoff,boxid);

%%%%%%%%%% Calculate movement onsets
if length(unique(boxid))<4, error('not bimanual task'); end % this isnt a bimanual task (bc most likely)

targeton = boxid==2; % left trials
lefttrials = u.findOnsetsAndOffsets(targeton); % these are the trials

targeton = boxid==3; % right trials
righttrials = u.findOnsetsAndOffsets(targeton); % these are the trials

accdata = [fdata(:,1), fdata(:,4)];  %left is 1, right is 4

tfcorrectonsets = double([0; diff(correctresponse(:)>1)>=1;]); % logical correct trial
wherecorrect = find(tfcorrectonsets)-1; % indeces where correct responses were dealt
righttrialsuccess = ismember(righttrials(:,2), wherecorrect); % right trials that are correct (where correct that dont align are experimentor dealt)
lefttrialsuccess = ismember(lefttrials(:,2), wherecorrect); % left trials that are correct
