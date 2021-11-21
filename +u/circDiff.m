function d = circDiff(data1,data2)
    d = pi - abs(pi - (mod(data1-data2,2*pi)));
end