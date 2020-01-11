function [ x_k,error ] = find_point_x_k(solution_NLP,solution_LP,Q,r,l)
%FIND_POINT_X_K 此处显示有关此函数的摘要
%   此处显示详细说明
% A_emission=reshape(repmat(c,1,T)',1,N*T)*power(delta_Pit,2);
% B_emission=0;
% for i=1:N
%     for t=1:T
%         B_emission=B_emission+2*c(i)*delta_Pit((i-1)*T+t,1)*Pit_LP((i-1)*T+t,1)+...
%             b(i)*delta_Pit((i-1)*T+t,1)+...
%             a(i)*delta_Uit((i-1)*T+t,1);
%     end
% end
% B_emission=B_emission-delta_emission_buy+delta_emission_sell;
% 
% C_emission=0;
% for i=1:N
%     for t=1:T
%         C_emission=C_emission+a(i)*Uit_LP((i-1)*T+t,1)+c(i)*Pit_LP((i-1)*T+t,1)^2+...
%             b(i)*Pit_LP((i-1)*T+t,1);
%     end
% end
% C_emission=C_emission-emission_buy_LP+emission_sell_LP-E0;

deta_x=solution_NLP-solution_LP;
a=deta_x'*Q*deta_x;
b=solution_LP'*Q*deta_x+deta_x'*Q*solution_LP+l'*deta_x;
c=solution_LP'*Q*solution_LP+l'*solution_LP-r;
[lambda]=solve_equation(a,b,c);

if(lambda~=-100)
    
    x_k=lambda*deta_x+solution_LP;
    error=0;
    max_g=x_k'*Q*x_k+l'*x_k-r;
    
 
else
    error=1;
    disp('λ错误');
end
end

