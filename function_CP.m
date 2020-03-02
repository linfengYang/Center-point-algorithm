function [ x,result,time ] = function_CP(  model,fid )
%UNTITLED 此处显示有关此函数的摘要
%model represents UC-CET formulation,fid represents the filepoint of log file
%miu=1/(1+1e3*exp(-3)^n)

%init
TimeLimit=model.TimeLimit;
tolerances.TimeLimit=TimeLimit;
solution_limit=1;
tolerances.solution_limit=solution_limit;
tolerances.heuristicfreq=0;
termination_LP=1e-3;
total_time=0;
iteration=1000;
termination_r=1e-3;

iter_num=0;
L=model.L;
location_uit=model.location_uit;

model_orgin.H=model.H;
model_orgin.f=model.f;
model_orgin.Q=model.Q;
model_orgin.l=model.l;
model_orgin.r=model.r;
model_orgin.col_H=size(model_orgin.H,1);

[ model ] = obj_linear( model,L );%Objective function linearization

[ new_model ] = extend_matric( model,1,-inf,inf,'C');%extend matrix,equivalent to add new variable

%Split CET constraints
Variable_num=length(new_model.ctype);
Aineq_CET=sparse(1,Variable_num);
Aineq_CET(1,Variable_num)=1;
Aineq_CET=Aineq_CET+new_model.l';
bineq_CET=new_model.r;
new_model.Aineq=[new_model.Aineq;Aineq_CET];
new_model.bineq=[new_model.bineq;bineq_CET];
two_bin_location=strfind(new_model.ctype,'B');
loacation_pit=find(ones(1,Variable_num)*new_model.Q~=0);

%NLP
ctype=model.ctype;
ctype(1,two_bin_location)='C';
tolerances.TOLER=5e-3;
[result,x,time ] = sovle( 2*model.Q,model.l,model.Aineq,model.bineq,model.Aeq,model.beq,[],[],[],model.lb,model.ub,ctype,5,tolerances);
x_k_1=x;
solution_NLP=x;% x_nlp
total_time=total_time+time;
TimeLimit=TimeLimit-time;
[max_g,obj]=inform(model_orgin,x);%Get the value of the objective function and CET constraints g(x)
disp_result( fid,'NLP',0,max_g,total_time,result,obj,0);
%Generate a good linear approximation for g(x)
[x,Aineq_PW,bineq_PW,iter_num,total_time,TimeLimit]=Step_Linear_approximation(new_model,model_orgin,fid,solution_NLP,total_time,TimeLimit,iter_num,termination_LP,iteration,Variable_num,location_uit,loacation_pit);
%Find_integer_solutions
[x,result,time,max_g]=Step_Find_integer_solutions(new_model,model_orgin,fid,Aineq_PW,bineq_PW,iter_num,solution_NLP,total_time,TimeLimit,location_uit,loacation_pit);

end

function [x,Aineq_PW,bineq_PW,iter_num,total_time,TimeLimit]=Step_Linear_approximation(model,model_orgin,fid,solution_NLP,total_time,TimeLimit,iter_num,termination_LP,iteration,Value_num,location_uit,loacation_pit)
%init
param.TOLER=5e-3;
param.solution_limit=0;
Aineq_PW=[];
bineq_PW=[];
x_k_1=[solution_NLP;0];
two_bin_location=strfind(model.ctype,'B');
ctype=model.ctype;
ctype(1,two_bin_location)='C';
max_g=1;
iter_num=1;

while(max_g>termination_LP&&iter_num<iteration)
    if(iter_num>1)
        %linear search
        [x_k,error]=find_point_x_k([solution_NLP;0],solution_LP,model.Q,model.r,model.l);
        if(error==0)
            %construct perspective-cut
            gradient=sparse(2*x_k'*model.Q);
            gradient(1,Value_num)=-1;
            uit=-(x_k.*x_k)'.*(ones(1,Value_num)*model.Q);
            gradient(1,location_uit)=uit(1,loacation_pit);
            
            bineq_perspective_cut=0;
            Aineq_perspective_cut=gradient;
            %add new perspective_cut
            Aineq_PW=[Aineq_PW;Aineq_perspective_cut];
            bineq_PW=[bineq_PW;bineq_perspective_cut];
            %solve LP problem
            param.TimeLimit=TimeLimit;
            [result,x,time ] = sovle(model.H,model.f,[model.Aineq;Aineq_PW],[model.bineq;bineq_PW],model.Aeq,model.beq,[],[],[],model.lb,model.ub,ctype,1,param);
            
        else
            disp('无法找到X_k!');
        end
    else
        param.TimeLimit=TimeLimit;
        [result,x,time ] = sovle( model.H,model.f,model.Aineq,model.bineq,model.Aeq,model.beq,[],[],[],model.lb,model.ub,ctype,1,param);
    end
    deta_x=((x-x_k_1)'*(x-x_k_1))^0.5;
    x_k_1=x;
    solution_LP=x;
    total_time=total_time+time;
    TimeLimit=TimeLimit-time;
    [max_g,obj]=inform(model,x);%get objective value and g(x)
    disp_result( fid,'LP',iter_num,max_g,total_time,result,obj,deta_x);
    iter_num=iter_num+1;
end
end


function [x,result_obj,time_find_solutions,g_x]=Step_Find_integer_solutions(model,model_orgin,fid,Aineq_PW,bineq_PW,iter_num,solution_NLP,total_time,TimeLimit,location_uit,loacation_pit)
%Judging whether the algorithm is timed out
if(TimeLimit<0)
    time_find_solutions=total_time;
    result_obj=null;
    return;
end
%init
Aineq_obj_cutting=[];
bineq_obj_cutting=[];
params.TOLER= 1e-3;
CC_iteration=1;
%extend model
[ new_model ] = extend_matric( model,1,0,inf,'C');
two_bin_location=strfind(model.ctype,'B');
ctype_NLP=model.ctype;
ctype_NLP(1,two_bin_location)='C';

Value_num=length(new_model.ctype);%get number of variables of new extend model
Value_num_NLP=length(ctype_NLP); %get number of variables of NLP problem model
%reformulate objective function
f_MILP=sparse(Value_num,1);
f_MILP(Value_num,1)=-1;

while(TimeLimit>0)
    if(CC_iteration>1)
        if(g_x<=0.001)
            %Fixed-integer neighborhood search step
            Aeq_UC_init_test=sparse(1:length(two_bin_location),two_bin_location,1,length(two_bin_location ),Value_num_NLP);
            beq_UC_init_test=x(two_bin_location,1);
            params.solution_limit=0;
            params.TOLER= 1e-3;
            params.TimeLimit=TimeLimit;
            [result,x,time ] = sovle( model.H,model.f,model.Aineq,model.bineq,[model.Aeq;Aeq_UC_init_test],[model.beq;beq_UC_init_test],model.Q,model.l,model.r,model.lb,model.ub,ctype_NLP,3,params);
            if(isempty(x))
                x=x_k_1;
                [g_x,obj]=inform(model,x);
                time_find_solutions=total_time;
                result_obj=obj;
                return;
            end
            
            deta_x=((x(1:model_orgin.col_H,1)-x_k_1(1:model_orgin.col_H,1))'*(x(1:model_orgin.col_H,1)-x_k_1(1:model_orgin.col_H,1)))^0.5;
            x_k_1=x;
            total_time=total_time+time;
            TimeLimit=TimeLimit-time;
            [g_x,obj]=inform(model,x);
            disp_result( fid,'neighborhood search NLP',iter_num,g_x,total_time,result,obj,deta_x);
            if(TimeLimit<0)
                time_find_solutions=total_time;
                result_obj=obj;
                return;
            end
            %X_Obj
            gradient=sparse(model.f');
            bineq_obj_cutting=gradient*x;
            Aineq_obj_cutting=gradient;

            if(g_x>-0.001&&g_x<0.001)
                gradient=sparse(2*x'*model.Q);
                uit=-(x.*x)'.*(ones(1,Value_num_NLP)*model.Q);
                gradient(1,location_uit)=uit(1,loacation_pit);
                gradient(1,Value_num_NLP)=-1;
                
                Aineq_PW=[Aineq_PW;gradient];
                bineq_PW=[bineq_PW;0];
            end

        else
            %Integer variable feasibility judgment step
            Aeq_UC_init_test=sparse(1:length(two_bin_location),two_bin_location,1,length(two_bin_location ),Value_num);
            beq_UC_init_test=x(two_bin_location,1);
            x_k_1=x;
            params.TimeLimit=TimeLimit;
            params.TOLER=1e-3;
            params.solution_limit=0;
            ctype=new_model.ctype;
            ctype(1,two_bin_location)='C';
            NLP_f(Value_num,1)=1;
            l=new_model.l;
            l(Value_num,1)=-1;
            [result,x,time ] = sovle( [],NLP_f,new_model.Aineq,new_model.bineq,...
                [new_model.Aeq;Aeq_UC_init_test],[new_model.beq;beq_UC_init_test],...
                new_model.Q,l,new_model.r,new_model.lb,new_model.ub,ctype,3,params);
            total_time=total_time+time;
            TimeLimit=TimeLimit-time;
            [g_x,obj]=inform(model,x);
            disp_result( fid,'NLP_feasible',iter_num,g_x,total_time,result,obj,0);
            if(result>1e-5)
                if(length(x_k_1)>length(solution_NLP))
                    x_k_1=x_k_1(1:length(solution_NLP)+1,1);
                else
                    x_k_1=[x_k_1;sparse(length(solution_NLP)-length(x_k_1),1)];
                end
                %linear search
                [x_k,error]=find_point_x_k([solution_NLP;0],x_k_1, model.Q, model.r, model.l);
                if(error==0)
                    gradient=sparse(2*x_k'*model.Q);
                    uit=-(x_k.*x_k)'.*(ones(1,length(model.l))*model.Q);
                    gradient(1,location_uit)=uit(1,loacation_pit);
                    gradient(1,length(model.l))=-1;
                    
                    Aineq_PW=[Aineq_PW;gradient];
                    bineq_PW=[bineq_PW;0];
                    
                else
                    disp('MILP x_k');
                end
            else
             %Fixed-integer neighborhood search step
                Aeq_UC_init_test=sparse(1:length(two_bin_location),two_bin_location,1,length(two_bin_location ),Value_num_NLP);
                beq_UC_init_test=x_k_1(two_bin_location,1);
                [result,x,time ] = sovle( [],model.f,model.Aineq,model.bineq,...
                    [model.Aeq;Aeq_UC_init_test],[model.beq;beq_UC_init_test],...
                    model.Q,model.l,model.r,model.lb,model.ub,ctype_NLP,3,params);
                
                if(isempty(x))
                    x=x_k_1;
                    [g_x,obj]=inform(model,x);
                    time_find_solutions=total_time;
                    result_obj=obj;
                    return;
                end
                
                deta_x=((x(1:model_orgin.col_H,1)-x_k_1(1:model_orgin.col_H,1))'*(x(1:model_orgin.col_H,1)-x_k_1(1:model_orgin.col_H,1)))^0.5;
                x_k_1=x;
                total_time=total_time+time;
                TimeLimit=TimeLimit-time;
                [g_x,obj]=inform(model,x);
                disp_result( fid,'neighborhood search NLP',iter_num,g_x,total_time,result,obj,deta_x);
                if(TimeLimit<0)
                    time_find_solutions=total_time;
                    result_obj=obj;
                    return;
                end
                
                
                gradient=sparse(model.f');
                bineq_obj_cutting=gradient*x;
                Aineq_obj_cutting=gradient;

                if(g_x>-0.001&&g_x<0.001)
                    gradient=sparse(2*x'*model.Q);
                    uit=-(x.*x)'.*(ones(1,Value_num_NLP)*model.Q);
                    gradient(1,location_uit)=uit(1,loacation_pit);
                    gradient(1,Value_num_NLP)=-1;
                    
                    Aineq_PW=[Aineq_PW;gradient];
                    bineq_PW=[bineq_PW;0];
                end
            end
        end
        %The integer ellipsoid center-point step
        params.solution_limit=1;
        params.TimeLimit=TimeLimit;
        params.TOLER=1e-3;
        miu=1/(1+1e3*exp(-3)^(CC_iteration-1));

        for i=1:length(bineq_PW)
            Aineq_PW_cutting(i,:)=[Aineq_PW(i,:),miu*(Aineq_PW(i,:)*Aineq_PW(i,:)')^0.5];
            bineq_PW_cutting(i,:)=bineq_PW(i,:);
        end
        for i=1:length(bineq_obj_cutting)
              Aineq_obj_cutting(1,Value_num)=0;
              Aineq_obj_cutting(1,Value_num)=(1/miu)*(Aineq_obj_cutting*Aineq_obj_cutting')^0.5;%
        end
        [result,x,time ] = sovle( [],f_MILP,[new_model.Aineq;Aineq_PW_cutting;Aineq_obj_cutting],[new_model.bineq;bineq_PW_cutting;bineq_obj_cutting],new_model.Aeq,new_model.beq,[],[],[],new_model.lb,new_model.ub,new_model.ctype,2,params);
        x_k_1=[x_k_1;0];
    else
        %first iteration
        params.solution_limit=2;
        params.TimeLimit=TimeLimit;
        miu=1/(1+1e3*exp(-3)^CC_iteration);
        for i=1:length(bineq_PW)
            Aineq_PW_cutting(i,:)=[Aineq_PW(i,:),miu*(Aineq_PW(i,:)*Aineq_PW(i,:)')^0.5];
            bineq_PW_cutting(i,:)=bineq_PW(i,:);
        end
        [result,x,time ] = sovle( [],f_MILP,[new_model.Aineq;Aineq_PW_cutting],[new_model.bineq;bineq_PW_cutting],new_model.Aeq,new_model.beq,[],[],[],new_model.lb,new_model.ub,new_model.ctype,2,params);
        x_k_1=x;
    end
    if(isempty(x))
        x=x_k_1;
        time_find_solutions=total_time;
        result_obj=obj;
        [g_x,obj]=inform(model,x);
        break;
    end
    total_time=total_time+time;
    deta_x=((x-x_k_1)'*(x-x_k_1))^0.5;
    x_k_1=x;
    [g_x,obj]=inform(model,x);
    TimeLimit=TimeLimit-time;
    disp_result( fid,'integer ellipsoid center-point MILP',iter_num,g_x,total_time,result,obj,deta_x);
    iter_num=iter_num+1;
    CC_iteration=CC_iteration+1;
end


end

function [max_g,obj]=inform(model,x)
    model.col_H=length(model.f);
    max_g=x(1:model.col_H,1)'*model.Q*x(1:model.col_H,1)+model.l'*x(1:model.col_H,1)-model.r;
    if(isempty( model.H))
        obj=model.f'*x(1:model.col_H,1);
    else
        obj=0.5*x(1:model.col_H,1)'*model.H*x(1:model.col_H,1)+model.f'*x(1:model.col_H,1);
    end
end