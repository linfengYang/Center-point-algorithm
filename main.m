%main
%Li Wei
%GuangXi university

clear all;

times_T=1;
times_Uit=1;
L =4;            %piecewise parameter
TimeLimit=3600;
FileFold='UC_AF\';

%Get filename of all instance
FileFoder=fullfile('.\UC_AF\');
dirOutput=dir(fullfile(FileFoder,'c*.mod'));
List_fileName={dirOutput.name}';
File_num=size(List_fileName,1); % total number of file
index=1;
while(index<=File_num)
    %Generate file path
    FileName=List_fileName{index};
    FileName=strcat(FileName,'');
    FilePath=strcat(FileFoder,FileName);
    
    dataUC=readdataUC(FilePath);      %Read file data
    
    FileType='.mod';
    fid_CP=fopen(strcat('result_',FileFold,FileName,'_CP_',num2str(times_Uit),'_times',FileType),'wt');
    %fid_CPLEX=fopen(strcat('result_',FileFold,FileName,'_CPLEX_',num2str(times_Uit),'_times',FileType),'wt');
    
    model=constraints_produce(dataUC); %generate UC_CET formulation
    
    H_num=size(model.H,1);
    H=model.H;
    f=model.f;
    Q=model.Q;
    l=model.l;
    r=model.r;
    
    %Center-point algorithm
    model.L=L;
    model.TimeLimit=TimeLimit;
    [ x,result,time ] = function_CP( model,fid_CP);
    
    %If you want to call cplex, remove the comment 
    %cplex
%     [ new_model_1 ] = CET_linear_piecewise_perspective_cut( model,4 );
%     [ model ] = obj_linear( model,4 );
%     cplex_solver = Cplex('Master problem for HTC');
%     cplex_solver.Model.sense = 'minimize';
%     cplex_solver.Model.Q = [];
%     cplex_solver.Model.obj = model.f;
%     cplex_solver.Model.lb = model.lb;
%     cplex_solver.Model.ub = model.ub;
%     cplex_solver.Model.A = [model.Aineq;model.Aeq];
%     cplex_solver.Model.lhs = [-Inf.*ones(size(model.bineq,1),1);model.beq];
%     cplex_solver.Model.rhs = [model.bineq;model.beq];
%     cplex_solver.Model.ctype = ctype;
%     cplex_solver.addQCs(model.l, model.Q,'L',model.r);
%     cplex_solver.InfoCallback.func = @mipex4cb;%调用callback
%     cplex_solver.InfoCallback.data.fid=fid_CPLEX_linear;
%     cplex_solver.InfoCallback.data.Q=model.Q;
%     cplex_solver.InfoCallback.data.r=model.r;
%     cplex_solver.InfoCallback.data.l=model.l;
%     cplex_solver.Param.timelimit.Cur=TimeLimit;
%     cplex_solver.DisplayFunc= 'off'; %不显示Cplex求解过程
%     cplex_solver.solve();
%     if(cplex_solver.Solution.status==108)
%         fprintf(fid_CPLEX_linear,'超时');
%     else
%         result_cplex=cplex_solver.Solution.objval;
%         gap_cplex=cplex_solver.Solution.miprelgap;
%         relaxa_solution=cplex_solver.Solution.bestobjval;
%         x=cplex_solver.Solution.x;
%         time=cplex_solver.Solution.time;
%         result_Obj_cplex=0.5*x(1:H_num,1)'*H*x(1:H_num,1)+f'*x(1:H_num);
%         CET_cplex=x(1:H_num,1)'*Q*x(1:H_num,1)+l'*x(1:H_num,1)-r;
%         fprintf(fid_CPLEX,'%d    %f    %f     %f\n',0,time,result_Obj_cplex,CET_cplex);
%     end

    index=index+1;
end

% function stop = mipex4cb(info,data)
% if info.IncObj<1e73
%     gap = info.MipGap * 100.0;
%     info.Time
%     info.BestObj
%     info.IncObj
%     if(~isempty(info.IncX))
%         x=info.IncX;
%         CET_cplex_linear=x'*data.Q*x+data.l'*x-data.r;
%         fprintf(data.fid,'%f    %f    %f     %f\n',info.Time,info.BestObj,info.IncObj,CET_cplex_linear);
%     end
% end
% stop = false;
% end