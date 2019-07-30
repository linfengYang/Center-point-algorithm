function [ new_model ] = obj_linear( model,L )
%OBJ_LINEAR 此处显示有关此函数的摘要
%   此处显示详细说明
variation_num=size(model.ctype,2);
secondary_vari_location=find(diag(model.H)~=0);
secondary_vari_num=size(secondary_vari_location,1);
new_model_variation_num=variation_num+secondary_vari_num;
Aineq=[];
bineq=[];
for i=1:secondary_vari_num
    for j=0:L
        gradiant=sparse(1,new_model_variation_num);
        gradiant(1,i)=model.f(i,1)-0.5*model.H(secondary_vari_location(i),secondary_vari_location(i))*(j/L)^2;
        gradiant(1,secondary_vari_location(i))=model.H(secondary_vari_location(i),secondary_vari_location(i))*(j/L)+model.f(secondary_vari_location(i),1);
        gradiant(1,variation_num+i)=-1;
        Aineq=[Aineq;gradiant];
        bineq=[bineq;0];
    end
end

new_model.H=[];
model.f(1:secondary_vari_num,1)=0;
model.f(secondary_vari_location,1)=0;
new_model.f=[model.f;ones(secondary_vari_num,1)];
new_model.Aineq=[model.Aineq,sparse(size(model.Aineq,1),secondary_vari_num)];
new_model.bineq=model.bineq;
new_model.Aineq=[new_model.Aineq;Aineq];
new_model.bineq=[new_model.bineq;bineq];
new_model.Aeq=[model.Aeq,sparse(size(model.Aeq,1),secondary_vari_num)];
new_model.beq=model.beq;
new_model.Q=[model.Q,sparse(size(model.Q,1),secondary_vari_num);sparse(secondary_vari_num,size(model.Q,2)+secondary_vari_num)];
new_model.l=[model.l;sparse(secondary_vari_num,1)];
new_model.r=model.r;
new_model.lb=[model.lb;-inf*ones(secondary_vari_num,1)];
new_model.ub=[model.ub;inf*ones(secondary_vari_num,1)];
add_type(1,1:secondary_vari_num)='C';
new_model.ctype=strcat(model.ctype,add_type);
end

