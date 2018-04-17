function settings = addExchangeMedum(settings)
%Minimal media simulation
inMedia = {
    'glycogen[s]'
    'palmitate[s]'
    'H2O[s]'
    'CO2[s]'
    'O2[s]'
};
constraints =  [
    -1000
    -1000
    -1000
    -1000
    -0
];

outMedia = {
    'H2O[s]'
    'palmitate[s]'
    'glycogen[s]'
    'O2[s]'
    'CO2[s]'
    'obectiveMetabolite[s]'
    'human_biomass[sm3]'
};


outConstraints = [
    1000
    1000
    1000
    1000
    0
    1000
    1000
    ];

settings.primaryObjective = 'MuscleATP';
settings.inMedia = inMedia;
settings.inValues = constraints;
settings.outMedia = outMedia;
settings.outConstraints = outConstraints;

end

