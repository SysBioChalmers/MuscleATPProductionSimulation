function result = molToW(data)
    %Converts mMol/h to ml/min
    %1 ?mol O2 = .022391 ml
    result = (0.022391 * 10^6)*data;
    result = result / 60;
end

