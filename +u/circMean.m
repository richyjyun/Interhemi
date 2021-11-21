function m = circMean(data1)
    m = angle(nansum(exp(1j*data1)));
end