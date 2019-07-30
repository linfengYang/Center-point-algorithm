function [ new_model ] = extend_matric( model,n,add_lb,add_ub,add_ctype)
%EXTEND_MATRIC 此处显示有关此函数的摘要
%   此处显示详细说明
if(isempty(model.H))
    new_model.H=model.H;
else
    new_model.H=[model.H,sparse(size(model.H,1),1);sparse(1,size(model.H,2)+n)];
end
new_model.f=[model.f;sparse(n,1)];
new_model.Aineq=[model.Aineq,sparse(size(model.Aineq,1),n)];
new_model.bineq=model.bineq;
new_model.Aeq=[model.Aeq,sparse(size(model.Aeq,1),n)];
new_model.beq=model.beq;
if(isempty(model.l))
    new_model.l=[];
else
    new_model.l=[model.l;sparse(n,1)];
end
if(isempty(model.Q))
    new_model.Q=model.Q;
else
    new_model.Q=[model.Q,sparse(size(model.Q,1),n);sparse(n,size(model.Q,2)+n)];
end
if(isempty(model.r))
     new_model.r=[];
else
     new_model.r=model.r;
end

 new_model.lb=[model.lb;add_lb*ones(n,1)];
 new_model.ub=[model.ub;add_ub*ones(n,1)];
 add_ctype_str(1,1:n)=add_ctype;
 new_model.ctype=strcat(model.ctype,add_ctype_str);
end

