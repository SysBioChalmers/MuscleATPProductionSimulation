function result = molToW(data)
% http://book.bionumbers.org/how-much-energy-is-released-in-atp-hydrolysis/
% Theoretical KJ/mol ATP in muscle is arround 70
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2342936/
% In practice one gets arround 27-50 KJ/ATP due to mechanical inefficency
%Method 2:
%Kcal/mol ATP 28.1250 -> 117.675 KJ/mol ATP
%Net Mechanical efficency ~0.30% -> 35.3025 from:
%Total power output generated during dynamic knee
%extensor exercise at different contraction frequencies
%http://link.springer.com/article/10.1007%2Fs00421-009-1008-7
%Efficiency in cycling: a review
%Gertjan Ettema Æ Ha°vard Wuttudal Lora°s

result = data/3600; %from h to s = W
%result = result*35.3; %convert to KJ/s
result = result*30; %convert to KJ/s



end

