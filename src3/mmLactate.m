function rate = mmLactate(substrate, vmax)
    km = 10.73;
    rate = (vmax*substrate)./(km + substrate);
end

