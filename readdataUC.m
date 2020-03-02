function dataUC=readdataUC(pathAndFilename)
fprintf('\n\n\n');
disp('------------------------读取机组组合问题数据文件----------------------------');
%pathAndFilename
%% 1.读取电力系统参数
% [filename,path] = uigetfile('.mod','请选择数据文件');
% pathAndFilename = strcat(path,filename);
% pathAndFilename='UC_AF/10_0_1_w.mod';
% pathAndFilename='UC_AF/10_std.mod';
[Pathstr,Filename,Ext] = fileparts(pathAndFilename);
dataUC.pathAndFilename = pathAndFilename;
fid=fopen(pathAndFilename);
%忽略第1行数据
fgetl(fid);                    %ProblemNum
%读取第2行数据
tmp=fscanf(fid,'%s',1);                     %HorizonLen
dataUC.T=fscanf(fid,'%d',1);                       %HorizonLen   即：时段T
%读取第3行数据
tmp=fscanf(fid,'%s',1);                     %NumThermal
dataUC.N=fscanf(fid,'%d',1);                       %NumThermal   即：机组数N
% 忽略第4-10行数据
% NumHydro , NumCascade , LoadCurve , MinSystemCapacity, MaxSystemCapacity,
% MaxThermalCapacity, Loads
fgetl(fid);  % 上一行的结尾
for i = 1:7 
    fgetl(fid); 
end
% 读取第11行数据 负荷
dataUC.PD=fscanf(fid,'%f',[1,dataUC.T]);    dataUC.PD = dataUC.PD';                   % 负荷
% 忽略第12行数据
fgetl(fid);fgetl(fid);
% 读取第13行数据 备用
dataUC.spin = fscanf(fid,'%f',[1,dataUC.T]);  dataUC.spin = dataUC.spin';     % 备用
% 忽略第14行数据
fgetl(fid);fgetl(fid);
%读取N台机组参数
for i = 1:dataUC.N
    if size(strfind(dataUC.pathAndFilename,'_std.mod') ,1) ~= 0
        unit_parameters(i,:) = fscanf(fid,'%f',[1,17]);     
    else
        unit_parameters(i,:) = fscanf(fid,'%f',[1,16]);
    end
    fscanf(fid,'%s',1);
    dataUC.p_rampup(i) = fscanf(fid,'%f',1);
    dataUC.p_rampdown(i) = fscanf(fid,'%f',1);
end
fgetl(fid);fgetl(fid);fgetl(fid);

dataUC.p_rampup = dataUC.p_rampup';   
dataUC.p_rampdown= dataUC.p_rampdown';

dataUC.alpha = unit_parameters(:,4);
dataUC.beta = unit_parameters(:,3);
dataUC.gamma = unit_parameters(:,2);
dataUC.p_low = unit_parameters(:,5);
dataUC.p_up = unit_parameters(:,6);
dataUC.time_on_off_ini = unit_parameters(:,7);
dataUC.time_min_on = unit_parameters(:,8);
dataUC.time_min_off = unit_parameters(:,9);
fixedCost4startup = unit_parameters(:,14);
dataUC.p_initial = unit_parameters(:,16);


%Read CET data
N=dataUC.N;
FileName=strcat('CET\CET_',Filename,'.mod');%,'_',num2str(N)
CET=load(FileName);
dataUC.a=CET(1:N,3);
dataUC.b=CET(1:N,2);
dataUC.c=CET(1:N,1);
dataUC.emmission_buy_max=CET(N+1,1);
dataUC.emmission_sell_max=CET(N+1,2);
dataUC.E0=CET(N+1,3);
dataUC.price_buy=CET(N+2,1);
dataUC.price_sell=CET(N+2,2);


if size(strfind(dataUC.pathAndFilename,'_std.mod') ,1) ~= 0
    dataUC.Cold_hour = unit_parameters(:,17);
end

dataUC.u0 = sparse(dataUC.N,1);
dataUC.u0(find(dataUC.p_initial>0)) = 1;

dataUC.p_startup = dataUC.p_low;
dataUC.p_shutdown = dataUC.p_low;

if size(strfind(dataUC.pathAndFilename,'_std.mod') ,1) ~= 0
    dataUC.Cold_cost = 2*fixedCost4startup;
    if size(strfind(dataUC.pathAndFilename,'5_std.mod') ,1) ~= 0
        dataUC.Cold_cost = 1*fixedCost4startup;
    end
else 
    dataUC.Cold_hour = sparse(dataUC.N,1);
    dataUC.Cold_cost = 1*fixedCost4startup;    
end
dataUC.Hot_cost = fixedCost4startup;

if size(strfind(dataUC.pathAndFilename,'8_std.mod') ,1) ~= 0
    dataUC.PD = sum(dataUC.p_up) * [0.71; 0.65; 0.62;  0.6;  0.58;  0.58;  0.6;  0.64;
                                    0.73; 0.8;  0.82;  0.83;   0.82; 0.8;  0.79;  0.79;
                                    0.83; 0.91; 0.9;  0.88;   0.85;  0.84;  0.79;  0.74];
    dataUC.spin = 0.03 * dataUC.PD;
end

dataUC.projected = 0;
dataUC.expended = 0;

fclose(fid);



