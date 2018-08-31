function EE = O2toenergy(vO2, vCO2)
    RQ = vCO2./vO2;
    RQ(RQ>1) = 1;
    RQ(RQ<0.7) = 0.7;
    
    vier = @(VO2,RQ) VO2 .* ((((RQ-0.71)./0.29).*21.4)+(((1-RQ)./0.29).*19.4))/60;
    EE = vier(vO2, RQ);
end

