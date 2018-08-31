function result = molToW(data)
    result = data/3600; %from h to s = W
    result = result*27; %convert to KJ/s
end

