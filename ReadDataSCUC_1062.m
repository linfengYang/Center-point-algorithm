%% ***************************************************************
%               Author��Lingfeng Yang,Jiangyao Luo,Yan Xu,Zhenrong Zhang,Zhaoyang Dong
%               Email��ylf@gxu.edu.cn,landiljy@163.com
%               Start date:2017.12.7
%               Finish date:2017.12.28
%               Function description:Read the SCUC data
%% ***************************************************************
%%%  ���ļ�ȡ����
function SCUC_data=ReadDataSCUC_1062(FileName,File)
filedata=fopen(FileName);

fgetl(filedata);     %���Ե�һ��ע��
%%%  ��ȡ��2������
baseparameters = fscanf(filedata,'%f',[1,7]);
SCUC_data.baseparameters.busN          = baseparameters(1);              %%%  The total number of buses
SCUC_data.baseparameters.branchN        = baseparameters(2);              %%%  The total number of branches
SCUC_data.baseparameters.balanceBus         = baseparameters(3);              %%%  Reference bus
SCUC_data.baseparameters.standardP        = baseparameters(4);              %%%  Reference power
SCUC_data.baseparameters.iterationMax     = baseparameters(5);              %%% It is not used in this program; �ڵ㷨����������,�˴��ò���
SCUC_data.baseparameters.centerParameter  = baseparameters(6);              %%% It is not used in this program; ���Ĳ���,�˴��ò���
SCUC_data.baseparameters.unitN          = baseparameters(7);              %%%  The total number of unit
standardP = SCUC_data.baseparameters.standardP;


fgetl(filedata);     %����ע��
fgetl(filedata);     %����ע��
%%%  ��ȡ��������
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    paraG(k,:) = fscanf(filedata,'%f',[1,16]);
    k=k+1;
    checkZero  = fscanf(filedata,'%d',1);
end
SCUC_data.units.bus_G      = paraG(:,1);    %Index of unit
SCUC_data.units.alpha      = paraG(:,2);    %A constant term in fuel cost function of thermal unit
SCUC_data.units.beta       = paraG(:,3);    %The coefficient of the first degree in fuel cost function of thermal unit
SCUC_data.units.gamma      = paraG(:,4);    %The coefficient of the second degree in fuel cost function of thermal unit
SCUC_data.units.PG_up      = paraG(:,5)/standardP;  %Maximum power output of unit i
SCUC_data.units.PG_low     = paraG(:,6)/standardP;  %Minimum power output of unit i
SCUC_data.units.QG_up      = paraG(:,7)/standardP;   %It is not used in this program;�˴��ò���
SCUC_data.units.QG_low     = paraG(:,8)/standardP;   %It is not used in this program;�˴��ò���
SCUC_data.units.U_ini      = paraG(:,9);            % It is not used in this program;�˴��ò���
SCUC_data.units.T_off      = paraG(:,10);           % It is not used in this program;�˴��ò���
SCUC_data.units.T_on      = paraG(:,11);
SCUC_data.units.ramp       = paraG(:,12)/standardP; % Ramp up/dowm limit of unit i
SCUC_data.units.start_cost = paraG(:,13);   %It is not used in this program;����������,�˴��ò���
SCUC_data.units.N          = k - 1;       %The total number of unit


% %%%%%%%%%%%%%%%%%%%%%    Ŀ�꺯��ϵ���س˻�׼ֵ     %%%%%%%%%%%%%%%%%%%%
SCUC_data.units.alpha = SCUC_data.units.alpha/SCUC_data.baseparameters.standardP;                                             %
SCUC_data.units.gamma = SCUC_data.units.gamma*SCUC_data.baseparameters.standardP;
%
% %%%%%%%%%%%%%%%%%%%%%    Ŀ�꺯��ϵ���س˻�׼ֵ     %%%%%%%%%%%%%%%%%%%%



%%%  ��ȡ�ڵ��ѹ���½�����
fgetl(filedata);     %����ע��
fgetl(filedata);     %����ע��
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    paraV(k,:) = fscanf(filedata,'%f',[1,2]);
    k=k+1;
    checkZero  = fscanf(filedata,'%d',1);
end
if k > 1
    SCUC_data.busV.v_up  = paraV(:,1);     %%%It is not used in this program; �ڵ��ѹ�Ͻ�,�˴��ò���
    SCUC_data.busV.v_low = paraV(:,2);     %%%It is not used in this program; �ڵ��ѹ�½�,�˴��ò���
else
    SCUC_data.busV.v_up  = [];     %%%It is not used in this program; �ڵ��ѹ�Ͻ�,�˴��ò���
    SCUC_data.busV.v_low = [];     %%%It is not used in this program; �ڵ��ѹ�½�,�˴��ò���
end

%%%  ��ȡ֧·���� (��������ѹ��֧· )
fgetl(filedata);     %����ע��
cs = fgetl(filedata);     %����ע��
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    branch(k,:) = fscanf(filedata,'%f',[1,7]); %�����е�ֵ��ʾ���˶�������ͬ��֧·������
    k           = k+1;
    checkZero   = fscanf(filedata,'%d',1);
end
if k > 1
    SCUC_data.branch.I = branch(:,1);     %%%  Start point in branch ij without transformer
    SCUC_data.branch.J = branch(:,2);     %%%  End point in branch ij without transformer
    SCUC_data.branch.R = branch(:,4);     %%% It is not used in this program; ֧·����,�˴��ò���
    SCUC_data.branch.X = branch(:,5);     %%%  Reactance of branch ij without transformer
    SCUC_data.branch.B = branch(:,6);     %%%  It is not used in this program;֧·�Եص���,�˴��ò���
    SCUC_data.branch.P = branch(:,7)/standardP;   %%%It is not used in this program; ֧·���书���Ͻ�,�˴��ò���
else
    SCUC_data.branch.I = [];
    SCUC_data.branch.J = [];
    SCUC_data.branch.R = [];
    SCUC_data.branch.X = [];
    SCUC_data.branch.B = [];
    SCUC_data.branch.P = [];
end

%%%   ��ȡ��ѹ��֧·����
fgetl(filedata);     %����ע��
fgetl(filedata);     %����ע��
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    branchTransformer(k,:) = fscanf(filedata,'%f',[1,8]);
    k                      = k+1;
    checkZero              = fscanf(filedata,'%d',1);
end
if k > 1
    SCUC_data.branchTransformer.I = branchTransformer(:,1);   %%%  Start point in branch ij with transformer
    SCUC_data.branchTransformer.J = branchTransformer(:,2);   %%%  End point in branch ij with transformer
    SCUC_data.branchTransformer.R = branchTransformer(:,4); %%% It is not used in this program;
    SCUC_data.branchTransformer.X = branchTransformer(:,5);    %%%  Reactance of branch ij with transformer
    SCUC_data.branchTransformer.K = 1./branchTransformer(:,8);  %
    SCUC_data.branchTransformer.P = branchTransformer(:,7)/standardP;  %%%It is not used in this program;
else
    SCUC_data.branchTransformer.I = [];
    SCUC_data.branchTransformer.J = [];
    SCUC_data.branchTransformer.R = [];
    SCUC_data.branchTransformer.X = [];
    SCUC_data.branchTransformer.K = [];
    SCUC_data.branchTransformer.P = [];
end


%%%  24ʱ�θ������ݼ���ת����
fgetl(filedata);     %����ע��
fgetl(filedata);     %����ע��
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    load(k,:) = fscanf(filedata,'%f',[1,7]);
    k         = k+1;
    checkZero = fscanf(filedata,'%d',1);
end
if k > 1
    SCUC_data.totalLoad.PD_T = load(:,1)/standardP;    %%%  System load demand
    SCUC_data.totalLoad.QD_T = load(:,2)/standardP;    %%It is not used in this program; �˴��ò���
    R_T  = load(:,5:7)/standardP;  %It is not used in this program;
    SCUC_data.totalLoad.R_T  = sum(R_T,2);                  %%It is not used in this program; �˴��ò���
    SCUC_data.totalLoad.T = size(SCUC_data.totalLoad.R_T,1); %The number of time period
else
    SCUC_data.totalLoad.PD_T = [];
    SCUC_data.totalLoad.QD_T = [];
    R_T  = [];
    SCUC_data.totalLoad.R_T  = [];
    SCUC_data.totalLoad.T = [];
end

%%%  �ڵ㸺������
fgetl(filedata);     %����ע��
fgetl(filedata);     %����ע��
checkZero=fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    nodeSequence(k) = checkZero;
    nodePQ(k,:)     = fscanf(filedata,'%f',[1,2]);
    k               = k+1;
    checkZero       = fscanf(filedata,'%d',1);
end

if k > 1
    SCUC_data.busLoad.bus_PDQR = nodeSequence';  %The node sequence
    SCUC_data.busLoad.node_P   = nodePQ(:,1)/standardP;  %The load in each node.
    SCUC_data.busLoad.node_Q   = nodePQ(:,2)/standardP;  %It is not used in this program;
    
    SCUC_data.busLoad.bus_P_factor = SCUC_data.busLoad.node_P/sum(SCUC_data.busLoad.node_P);  %%%  �ڵ㸺������   �����ڵ�ռ���ɵı�������������������ܸ����з��为�ɡ�
    SCUC_data.busLoad.bus_Q_factor = SCUC_data.busLoad.node_Q/sum(SCUC_data.busLoad.node_Q);  %It is not used in this program;
    
    SCUC_data.busLoad.node_P = SCUC_data.totalLoad.PD_T * SCUC_data.busLoad.bus_P_factor';
    SCUC_data.busLoad.node_Q = SCUC_data.totalLoad.QD_T * SCUC_data.busLoad.bus_Q_factor';  %It is not used in this program;
else
    SCUC_data.busLoad.bus_PDQR = [];
    SCUC_data.busLoad.node_P   = [];
    SCUC_data.busLoad.node_Q   = [];  %It is not used in this program;
    
    SCUC_data.busLoad.bus_P_factor = [];  %%%  �ڵ㸺������   �����ڵ�ռ���ɵı�������������������ܸ����з��为�ɡ�
    SCUC_data.busLoad.bus_Q_factor = [];  %It is not used in this program;
    
    SCUC_data.busLoad.node_P = [];
    SCUC_data.busLoad.node_Q = [];  %It is not used in this program;
end


N=SCUC_data.units.N;
FileName=strcat('CET\CET_',File);%,'_',num2str(N)
fid_CET=fopen(FileName);
CET= fscanf(fid_CET,'%f',[3,N]);
SCUC_data.a=CET(3,1:N)';
SCUC_data.b=CET(2,1:N)';
SCUC_data.c=CET(1,1:N)';
CET= fscanf(fid_CET,'%f',[1,5]);
SCUC_data.price_buy=CET(1);
SCUC_data.price_sell=CET(2);
SCUC_data.emmission_buy_max=CET(3);
SCUC_data.emmission_sell_max=CET(4);
SCUC_data.E0=CET(5);

%����������е����ñ���
clearvars -except SCUC_data FileName;
end


