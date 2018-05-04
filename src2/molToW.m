function result = molToW(data)
    result = data/3600; %from h to s = W
    result = result*26.9; %convert to KJ/s
end

