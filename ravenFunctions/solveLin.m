function solution = solveLin(model, silent)
    if nargin < 2
        silent = false;
    end

    % Setup the problem to feed to Matlab lin solver.
    prob=[];
    prob.c = model.c*-1;
    prob.a = vertcat(model.S, -model.S);
    
    if size(model.b,2)==1
        prob.b = vertcat(model.b(:,1), model.b(:,1));
    else
        prob.b = vertcat(model.b(:,2), model.b(:,1));
    end
    
    prob.blx = model.lb;
    prob.bux = model.ub;
    
    options = optimoptions('linProg', 'Display', 'off');
    
    [x,fval,exitflag,output] = linprog(prob.c, prob.a, prob.b, [], [], prob.blx, prob.bux, [], options);

    if exitflag == 1
        solution.x = x;
        solution.f = fval;
    else
        solution.x = 0;
        solution.f = 0;       
    end

    
end