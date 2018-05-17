function rate = mmLactate(substrate)
    vmax = 9;
    km = 10.73;
    rate = (vmax*substrate)./(km + substrate);
end

