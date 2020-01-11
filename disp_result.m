function disp_result( fid,type,iter_num,max_g,total_time,result,obj,deta_x)
%DISP_RESULT 此处显示有关此函数的摘要
%   此处显示详细说明
if(iter_num==1)
    disp('迭代次数=');
    disp(iter_num);
    disp('时间=');
    disp(total_time);
    disp('结果=');
    disp(result);
    disp('max_g');
    disp(max_g);
    fprintf(fid,'迭代次数     类型      时间    结果       g              obj                 deta_x\n');
    fprintf(fid,'  %d     %s      %f       %f          %f                  %f                %f\n',iter_num,type,total_time,result,max_g,obj,deta_x);
else
    disp('迭代次数=');
    disp(iter_num);
    disp('时间=');
    disp(total_time);
    disp('结果=');
    disp(result);
    disp('max_g');
    disp(max_g);
    disp('obj');
    disp(obj);
       disp('deta_x');
    disp(deta_x);
    fprintf(fid,'  %d     %s      %f       %f       %f       %f                %f\n',iter_num,type,total_time,result,max_g,obj,deta_x);
end
end

