%main
%Li Wei
%GuangXi university

clear all;

times_T=1;
times_Uit=1;
L =4;            %piecewise parameter
TimeLimit=3600;
for i=1:2

    if(i==1)
        FileFold='UC_AF\';
        %Get filename of all instance
        FileFoder=fullfile('.\UC_AF\');
        dirOutput=dir(fullfile(FileFoder,'c*.mod'));
    else
        FileFold='DC\';
        %Get filename of all instance
        FileFoder=fullfile('.\DC\');
        dirOutput=dir(fullfile(FileFoder,'SCUC*.txt'));
    end
    
    List_fileName={dirOutput.name}';
    File_num=size(List_fileName,1); % total number of file
    index=1;
    while(index<=File_num)
        %Generate file path
        FileName=List_fileName{index};
        FileName=strcat(FileName,'');
        FilePath=strcat(FileFoder,FileName);
        FileType='.mod';
        if(i==1)
            dataUC=readdataUC(FilePath);      %Read file data
            fid_GA=fopen(strcat('result\','result_','_GA_',FileName,num2str(times_Uit),'_times',FileType),'wt');
            fid_PSO=fopen(strcat('result\','result_','_PSO_',FileName,num2str(times_Uit),'_times',FileType),'wt');
            model=constraints_produce(dataUC);
        else
            if(index==1)
                dataUC=ReadDataSCUC_1062(FilePath,FileName);
            else
                dataUC=ReadDataSCUC(FilePath,FileName);
            end
            fid_CP=fopen(strcat('result\','result_','_CP_',FileName,num2str(times_Uit),'_times',FileType),'wt');
            fid_CPLEX=fopen(strcat('result\','result_','_CPLEX_',FileName,num2str(times_Uit),'_times',FileType),'wt');
            model=constraints_produce_DC(dataUC);
        end
        
        H_num=size(model.H,1);
        H=model.H;
        f=model.f;
        Q=model.Q;
        l=model.l;
        r=model.r;
        
        if(i==1)
            
            [ model ] = obj_linear( model,4 );
            [ new_model ] = extend_matric( model,1,-inf,inf,'C');
            Variable_num=length(new_model.ctype);
            Aineq_CET=sparse(1,Variable_num);
            Aineq_CET(1,Variable_num)=1;
            Aineq_CET=Aineq_CET+new_model.l';
            bineq_CET=new_model.r;
            new_model.Aineq=[new_model.Aineq;Aineq_CET];
            new_model.bineq=[new_model.bineq;bineq_CET];
            new_l(Variable_num,1)=-1;
            new_model.l=new_l;
            f=@(x)Ftc(x,new_model.f);
            c=@(x)CET(x,new_model.Q,new_model.l);
            A=[new_model.Aineq;new_model.Aeq;-1*new_model.Aeq];
            b=[new_model.bineq;new_model.beq;-1*new_model.beq];
            LB=new_model.lb;
            UB=new_model.ub;
            nvars=length(new_model.ctype);
            IntCon=find(new_model.ctype=='B');
            %GA
            options = optimoptions('ga','UseVectorized',true,'MaxTime',TimeLimit,'Display','iter','FunctionTolerance',1e-9,'UseParallel',true,'MaxStallGenerations',100);
            clc;
            diary(strcat('result\','log_','_GA_',FileName,num2str(times_Uit),'_times',FileType));
            diary on;
            [x,fval,exitflag,output]=ga(f,nvars,A,b,[],[],LB,UB,c,IntCon,options);
            if(isempty(fval)==0)
                fprintf(fid_GA,' %f\n',fval(1)+0);
            else
                fprintf(fid_GA,'error');
            end
            diary off;

            %PSO
            clc;
            diary(strcat('result\','log_','_PSO_',FileName,num2str(times_Uit),'_times',FileType));
            diary on;
            [x, fval,time] = function_PSO(new_model,dataUC,TimeLimit);
            fprintf(fid_PSO,' %f    %f\n',fval,time);
            diary off;
        else
            %Center-point algorithm
            model.L=L;
            model.TimeLimit=TimeLimit;
            [ x,result,time ] = function_CP( model,fid_CP);
            
            
            %If you want to call cplex, remove the comment
            %cplex
            %     [ new_model_1 ] = CET_linear_piecewise_perspective_cut( model,4 );
            [ model ] = obj_linear( model,4 );
            cplex_solver = Cplex('Master problem for HTC');
            cplex_solver.Model.sense = 'minimize';
            cplex_solver.Model.Q = [];
            cplex_solver.Model.obj = model.f;
            cplex_solver.Model.lb = model.lb;
            cplex_solver.Model.ub = model.ub;
            cplex_solver.Model.A = [model.Aineq;model.Aeq];
            cplex_solver.Model.lhs = [-Inf.*ones(size(model.bineq,1),1);model.beq];
            cplex_solver.Model.rhs = [model.bineq;model.beq];
            cplex_solver.Model.ctype = model.ctype;
            cplex_solver.Param.timelimit.Cur=TimeLimit;
            cplex_solver.addQCs(model.l, model.Q,'L',model.r);
            cplex_solver.InfoCallback.func = @mipex4cb;%调用callback
            cplex_solver.InfoCallback.data.fid=fid_CPLEX;
            cplex_solver.InfoCallback.data.Q=model.Q;
            cplex_solver.InfoCallback.data.r=model.r;
            cplex_solver.InfoCallback.data.l=model.l;
            cplex_solver.DisplayFunc= 'off'; %不显示Cplex求解过程
            cplex_solver.solve();
            if(cplex_solver.Solution.status==108)
                fprintf(fid_CPLEX,'超时');
            else
                result_cplex=cplex_solver.Solution.objval;
                gap_cplex=cplex_solver.Solution.miprelgap;
                relaxa_solution=cplex_solver.Solution.bestobjval;
                x=cplex_solver.Solution.x;
                time=cplex_solver.Solution.time;
                result_Obj_cplex=0.5*x(1:H_num,1)'*H*x(1:H_num,1)+f'*x(1:H_num);
                CET_cplex=x(1:H_num,1)'*Q*x(1:H_num,1)+l'*x(1:H_num,1)-r;
                fprintf(fid_CPLEX,'%d    %f    %f     %f\n',0,time,result_Obj_cplex,CET_cplex);
            end
        end
        index=index+1;
    end
end


function stop = mipex4cb(info,data)
if info.IncObj<1e73
    gap = info.MipGap * 100.0;
    info.Time
    info.BestObj
    info.IncObj
    if(~isempty(info.IncX))
        x=info.IncX;
        CET_cplex_linear=x'*data.Q*x+data.l'*x-data.r;
        fprintf(data.fid,'%f    %f    %f     %f\n',info.Time,info.BestObj,info.IncObj,CET_cplex_linear);
    end
end
stop = false;
end

function y=Ftc(x,f)
value_num=length(f);
y=x(:,1)*f(1);
for i=2:value_num
    y=y+x(:,i)*f(i);
end
end

function [c,ceq]=CET(x,Q,l)
value_num=length(l);
second_param=diag(Q);
c=x(:,1).^2*second_param(1)+x(:,1)*l(1);
for i=2:value_num
    c=c+x(:,i).^2*second_param(i)+x(:,i)*l(i);
end
ceq=[];
end

function [c]=CET_1(x,Q,l)
value_num=length(l);
second_param=diag(Q);
c=x(:,1).^2*second_param(1)+x(:,1)*l(1);
for i=2:value_num
    c=c+x(:,i).^2*second_param(i)+x(:,i)*l(i);
end

end



