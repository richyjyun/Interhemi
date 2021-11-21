function [data, fs, chnm, itrigch] = LoadGug(fname, dwn)
% Loads a guger bin file with relevant meta data from the cfg file.
% Based on a combination of parts of trainalign and loadgug_baseline
% fname is the file name with no extensions, dwn is the downsampling rate

if(nargin<2)
    dwn = 1;
end

load([fname '.cfg'],'-mat');

nchans = length(find(UI.ch_enabled))+length(find(UI.ga_trigger)); % number of channels
if isfield(UI, 'wd_channel')
    nchans = nchans + 1; % Window discriminator indicates daqdiscrim file.
else
    nchans = nchans + 2; % Assume daqbinmanual was used.
end
chnm = cell(nchans,1); % name of channels
chgu = zeros(nchans,1); % gUSBamp that channel was recorded on
ii = 0;
for iga = 1:length(UI.ga_trigger)
    ind = find(UI.ch_enabled(iga,:));
    chnm(ii+1:ii+length(ind)) = UI.ch_name(iga,ind)';
    chgu(ii+1:ii+length(ind)) = iga*ones(length(ind),1);
    ii = ii + length(ind);
    if UI.ga_trigger(iga), chnm{ii+1} = ['trig ' num2str(iga)]; chgu(ii+1) = iga; ii = ii + 1; end
end
if isfield(UI, 'wd_channel')
    chnm{end} = 'Discrim';  % One extra channel in daqdiscrim files.
else
    chnm{end-1} = 'Behave1'; % Two extra chanels in daqbimanual files.
    chnm{end} = 'Behave2';
end
itrigch = [strmatch('trig 1',chnm), strmatch('trig 2',chnm)];
fs = UI.samprate/dwn;

skip = nchans*(dwn-1)*4;
precision = sprintf('%d*single',nchans);
fid = fopen([fname,'.bin'], 'r');
fseek(fid,0,'bof');
data = fread(fid,[1,inf],precision,skip);%fread(fid,[nchans, inf],'*single')';
data = reshape(data,[nchans,length(data)/nchans])';
fclose(fid);

end