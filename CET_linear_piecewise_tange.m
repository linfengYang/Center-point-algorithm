function [ new_model ] = CET_linear_piecewise_tange(  model,L )
%CET_LINEAR_1 此处显示有关此函数的摘要
%   没有使用Zit而是一个变量Z，切线的表达式为一阶泰勒展开
variation_num=size(model.ctype,2);
secondary_vari_location=find(diag(model.H)~=0);
c=model.Q*ones(variation_num,1);
c=c(secondary_vari_location);
secondary_vari_num=size(secondary_vari_location,1);
new_model_variation_num=variation_num+1;

two_bin_vari_location=find(model.ctype=='B');
two_bin_vari_coefficence_l=model.l(two_bin_vari_location);

Aineq=[];
bineq=[];
for i=0:L
        gradiant=[sparse(1,secondary_vari_num),2*c'*(i/L),sparse(1,variation_num-2*secondary_vari_num),-1];
        Aineq=[Aineq;gradiant];
        bineq=[bineq;sum(c*(i/L)^2)];
end

Aineq_Z=model.l;
Aineq_Z(new_model_variation_num,1)=1;
bineq_Z=model.r;
Aineq=[Aineq;Aineq_Z'];
bineq=[bineq;bineq_Z];


new_model.H=[model.H,sparse(size(model.H,1),1);sparse(1,new_model_variation_num)];
new_model.f=[model.f;0];
new_model.Aineq=[model.Aineq,sparse(size(model.Aineq,1),1)];
new_model.bineq=model.bineq;
new_model.Aineq=[new_model.Aineq;Aineq];
new_model.bineq=[new_model.bineq;bineq];
new_model.Aeq=[model.Aeq,sparse(size(model.Aeq,1),1)];
new_model.beq=model.beq;
new_model.Q=[];
new_model.l=[];
new_model.r=[];
new_model.lb=[model.lb;-inf];
new_model.ub=[model.ub;inf];
add_type(1,1)='C';
new_model.ctype=strcat(model.ctype,add_type);

end