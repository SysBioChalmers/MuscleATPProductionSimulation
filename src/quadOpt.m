function solution = quadOpt(model)
    nrOfAttempts = 3;
    options.solver = 2;  
    options.maxnumiterations  = 8000;
    options.maxnumseconds = 1.5*60;
    options.primaltolerance  = 10^-5;
    options.dualtolerance    = 10^-5;
    options.verbose = 0;

    q = zeros(length(model.lb));
    objective = find(model.qc);
    objectiveVal = abs(model.qc(objective));
    for i = 1:length(objective)
        q(objective(i), objective(i)) = 2*objectiveVal(i);
    end
    c = -model.c;
    S = model.S;
    b = model.b;
    lb = model.lb;
    ub = model.ub;

    for i=1:nrOfAttempts
        [x,fval,flag] = clp(q,c,[],[],S,b,lb,ub,options);    
        if flag == 0
            break
        else
            disp(flag)
        end
    end
    
    if flag ~= 0
        solution.x = 0;
        solution.f = 0;
    else
        solution.x = x;
        solution.f = fval;
    end
end