function p = getPower(data,bp,fs)
p = abs(hilbert(bpfilt(data,bp,fs,3)));
p = mean(p,2);
end