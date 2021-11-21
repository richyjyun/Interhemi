function data = FilterAccForMovementTimes(data, fs, band)
% this function is intended to take accelerometer data and extract time points
% for which a movement has occured
%
% arb 10 nov 2014

if ~exist('band'), band = [2 20]; end 

if isempty(fs), data = []; return; end

if strcmp(band, 'richardson')
    bpf = [10 150]; % [10,50] chL/R transform: 1. bandpass filter (Hz)
    rms = 5; % 5 chL/R transform: 2. lowpass filter (Hz)
    [bbpf,abpf] = butter(1,bpf/(fs/2)); % 1st order bandpass (see daq_sapi_*)
    [brms,arms] = butter(1,rms/(fs/2),'low'); % 2nd order lowpass (see daq_sapi_*)
    data = filtfilt(bbpf,abpf,data(:)'); %'/1e3); % filter with bandpass
%     data = abs(data);
    data = sqrt(max(filtfilt(brms,arms,data.^2),zeros(1,length(data)))); % filter with lowpass
    return
end

Wn_theta = [band(1)/(fs/2) band(2)/(fs/2)]; % normalized by the nyquist frequency

[btheta,atheta] = butter(3,Wn_theta);

data = filtfilt(btheta,atheta,data);

data = abs(hilbert(data));

%     [bbpf,abpf] = butter(1,bpf/(UI.samprate/2)); % 1st order bandpass (see daq_sapi_*)
%     [brms,arms] = butter(2,rms/(UI.samprate/2),'low'); % 2nd order lowpass (see daq_sapi_*)
%      left = sqrt(max(filter(brms,arms,filter(bbpf,abpf,dat(iL,:)/1e3).^2),zeros(1,size(dat,2)))); % RMS accel left
%         right = sqrt(max(filter(brms,arms,filter(bbpf,abpf,dat(iR,:)/1e3).^2),zeros(1,size(dat,2)))); % RMS accel right
%         

%original richardson
% if strcmp(band, 'richardson')
%     bpf = [10 50]; % chL/R transform: 1. bandpass filter (Hz)
%     rms = 5; % chL/R transform: 2. lowpass filter (Hz)
%     [bbpf,abpf] = butter(1,bpf/(fs/2)); % 1st order bandpass (see daq_sapi_*)
%     [brms,arms] = butter(2,rms/(fs/2),'low'); % 2nd order lowpass (see daq_sapi_*)
%     data = filter(bbpf,abpf,data(:)'); %'/1e3); % filter with bandpass
%     data = sqrt(max(filter(brms,arms,data.^2),zeros(1,length(data)))); % filter with lowpass
%     return
% end