function [result,x,time ] = sovle( H,f,Aineq,bineq,Aeq,beq,Q,l,r,lb,ub,ctype,choice_solver,tolerances)
%  According to the value of choice_solver to choose different solve

toler=tolerances.TOLER;
TimeLimit=tolerances.TimeLimit;
solution_limit=tolerances.solution_limit;
% heuristicfreq=tolerances.heuristicfreq;
option = cplexoptimset;
option.Display = 'off';
option.mip.tolerances.mipgap=toler;
option.timelimit=TimeLimit;

if(solution_limit==0)
    %"solution_limit==0" means the number of solution limit in cplex is defaut
else
    option.mip.limits.solutions=solution_limit;
end

%Cplex class
cplex_class = Cplex('Milp for HTC');
cplex_class.Model.sense = 'minimize';
cplex_class.Model.Q = H;
cplex_class.Model.obj = f;
cplex_class.Model.lb = lb;
cplex_class.Model.ub = ub;
cplex_class.Model.A = [Aineq;Aeq];
cplex_class.Model.lhs = [-Inf.*ones(size(bineq,1),1);beq];
cplex_class.Model.rhs = [bineq;beq];
cplex_class.Model.ctype = ctype;
cplex_class.Param.timelimit.Cur=TimeLimit;
cplex_class.DisplayFunc= 'off';



switch choice_solver
    case 1
        cplex_class.solve();
        x = cplex_class.Solution.x;
        time=cplex_class.Solution.time;
        result=cplex_class.Solution.objval;
        
    case 2
        [x, result, exitflag, output] = cplexmilp (f, Aineq, bineq, Aeq, beq,...
            [ ], [ ], [ ], lb, ub, ctype, [ ], option);
        time=output.time;
        
    case 3
        [x, result, exitflag, output] = cplexqcp (H, f, Aineq, bineq, Aeq, beq, l, Q, r, lb, ub, ctype, option);
        time=output.time;
        
    case 4
        [x,result,exitflag,output]=cplexmiqcp(H, f, Aineq, bineq, Aeq, beq, l, Q, r, [], [], [], lb, ub, ctype, [], option);
        time=output.time;
    case 5
        cplex_class.Param.mip.tolerances.mipgap.Cur = toler;
        cplex_class.solve();
        x = cplex_class.Solution.x;
        time=cplex_class.Solution.time;
        result=cplex_class.Solution.objval;
end


disp('Time=');
disp(time);
fprintf ('Solution value = %f \n', result);

end

