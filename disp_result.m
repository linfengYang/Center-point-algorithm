function disp_result( fid,type,iter_num,max_g,total_time,result,obj,deta_x)
%DISP_RESULT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
if(iter_num==1)
    disp('��������=');
    disp(iter_num);
    disp('ʱ��=');
    disp(total_time);
    disp('���=');
    disp(result);
    disp('max_g');
    disp(max_g);
    fprintf(fid,'��������     ����      ʱ��    ���       g              obj                 deta_x\n');
    fprintf(fid,'  %d     %s      %f       %f          %f                  %f                %f\n',iter_num,type,total_time,result,max_g,obj,deta_x);
else
    disp('��������=');
    disp(iter_num);
    disp('ʱ��=');
    disp(total_time);
    disp('���=');
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

