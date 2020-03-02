%% ***************************************************************
%             ���ߣ�ylf
%             ԭ�����ڣ�2017��12��7��
%             �޸����ڣ�
%             ����˵������SCUC�����ļ���ȡ���������
%% ***************************************************************
%%%  ���ļ�ȡ����
function SCUCdata=ReadDataSCUC(FileName,File)
%% ����˵����******************************************************
% �������˵����FileName�������ļ���
% �������˵����SCUC_data.nodeNum:�ڵ�����  SCUC_data.branchNum:֧·��  standardP����׼����  iterationMax�����������
%centerParameter�����Ĳ���  IterationPrecision���������� objectFunc��Ŀ�꺯������ branchroute ��·�� BalanceI��ƽ��ڵ��  branchI����·ʼ�ڵ�i  branchJ����·�սڵ�j
% branchR����·����  branchX����·�翹  branchB����·�Եص���  branchGroundI���ӵ�֧·�ڵ�i  branchGroundB���ӵ�֧·����
% branchTransformerI:��ѹ���ڵ�i  branchTransformerJ����ѹ���ڵ�j  branchTransformerR����ѹ������  branchTransformerX����ѹ���翹
% branchTransformerK����ѹ�����  NodeI�����в����ڵ�i   NodePg:�й�����  NodeQg���޹�����  NodePl���й�����
% NodeQl:�޹�����  PvI��PV�ڵ�i  PvV:PV�ڵ��ѹ  PvQmin���޹���������  PvQmax���޹���������
% GenI:������ڵ��  GenC:�������������ϵ��  GenB:��������һ����ϵ��  GenA���������Զ�����ϵ��
% GenPmin���й���������  GenPmax���й���������  PvNum:Pv�ڵ���  GenNum��������ڵ���
%% ***************************************************************

filedata=fopen(FileName);

%fgetl(filedata);     %���Ե�һ��ע��
%%%  ��ȡ��2������
baseparameters = fscanf(filedata,'%f',[1,7]);
SCUCdata.baseparameters.busN          = baseparameters(1);              %%%  �ڵ�����
SCUCdata.baseparameters.branchN        = baseparameters(2);              %%%  ֧·��
SCUCdata.baseparameters.balanceBus         = baseparameters(3);              %%%  ƽ��ڵ�
SCUCdata.baseparameters.standardP        = baseparameters(4);              %%%  ��׼����
SCUCdata.baseparameters.iterationMax     = baseparameters(5);              %%%  �ڵ㷨����������
SCUCdata.baseparameters.centerParameter  = baseparameters(6);              %%%  ���Ĳ���
SCUCdata.baseparameters.unitN          = baseparameters(7);              %%%  ������
standardP = SCUCdata.baseparameters.standardP;


% fgetl(filedata);     %����ע��
% fgetl(filedata);     %����ע��
%%%  ��ȡ��������
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    paraG(k,:) = fscanf(filedata,'%f',[1,13]);
    k=k+1;
    checkZero  = fscanf(filedata,'%d',1);
end
SCUCdata.units.bus_G      = paraG(:,1);
SCUCdata.units.alpha      = paraG(:,2);
SCUCdata.units.beta       = paraG(:,3);
SCUCdata.units.gamma      = paraG(:,4);
SCUCdata.units.PG_up      = paraG(:,5)/standardP;
SCUCdata.units.PG_low     = paraG(:,6)/standardP;
SCUCdata.units.QG_up      = paraG(:,7)/standardP;
SCUCdata.units.QG_low     = paraG(:,8)/standardP;

SCUCdata.units.T_off      = paraG(:,9);           % ��Сͣ��ʱ��
SCUCdata.units.T_on      = paraG(:,10);
SCUCdata.units.U_ini      = paraG(:,11);            % ��ʼ��������ʱ�䣿�Ҳ���Ҫ
SCUCdata.units.ramp       = paraG(:,12)/standardP; % ����Լ��
SCUCdata.units.start_cost = paraG(:,13);
% SCUCdata.units.a      = paraG(:,16);
% SCUCdata.units.b       = paraG(:,15);
% SCUCdata.units.c      = paraG(:,14);
SCUCdata.units.N          = k - 1;                 % ������Ŀ


% %%%%%%%%%%%%%%%%%%%%%    Ŀ�꺯��ϵ���س˻�׼ֵ     %%%%%%%%%%%%%%%%%%%%
% beta  = beta*standardP;                                              %
% gamma = gamma*(standardP*standardP);                                 %
SCUCdata.units.alpha = SCUCdata.units.alpha/SCUCdata.baseparameters.standardP;                                             %
SCUCdata.units.gamma = SCUCdata.units.gamma*SCUCdata.baseparameters.standardP;
% SCUCdata.units.a = SCUCdata.units.a/SCUCdata.baseparameters.standardP;                                             %
% SCUCdata.units.c = SCUCdata.units.c*SCUCdata.baseparameters.standardP;
%
% %%%%%%%%%%%%%%%%%%%%%    Ŀ�꺯��ϵ���س˻�׼ֵ     %%%%%%%%%%%%%%%%%%%%



%%%  ��ȡ�ڵ��ѹ���½�����
% fgetl(filedata);     %����ע��
fgetl(filedata);     %����ע��
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    paraV(k,:) = fscanf(filedata,'%f',[1,2]);
    k=k+1;
    checkZero  = fscanf(filedata,'%d',1);
end
SCUCdata.busV.v_up  = paraV(:,1);     %%% �ڵ��ѹ�Ͻ�
SCUCdata.busV.v_low = paraV(:,2);     %%% �ڵ��ѹ�½�

%%%  ��ȡ֧·���� (��������ѹ��֧· )
% fgetl(filedata);     %����ע��
fgetl(filedata);     %����ע��
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    branch(k,:) = fscanf(filedata,'%f',[1,7]);
    k           = k+1;
    checkZero   = fscanf(filedata,'%d',1);
end
SCUCdata.branch.I = branch(:,1);     %%%  ֧·ʼ��
SCUCdata.branch.J = branch(:,2);     %%%  ֧·�յ�
SCUCdata.branch.R = branch(:,4);     %%%  ֧·����
SCUCdata.branch.X = branch(:,5);     %%%  ֧·�翹
SCUCdata.branch.B = branch(:,6);     %%%  ֧·�Եص���  ��Ӧ�����pi�͵�ֵ��·�ĶԵص���
SCUCdata.branch.P = branch(:,7)/standardP;   %%% ֧·���书���Ͻ�

%%%   ��ȡ��ѹ��֧·����
% fgetl(filedata);     %����ע��
fgetl(filedata);     %����ע��
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    branchTransformer(k,:) = fscanf(filedata,'%f',[1,8]);
    k                      = k+1;
    checkZero              = fscanf(filedata,'%d',1);
end
SCUCdata.branchTransformer.I = branchTransformer(:,1);
SCUCdata.branchTransformer.J = branchTransformer(:,2);
SCUCdata.branchTransformer.R = branchTransformer(:,4);
SCUCdata.branchTransformer.X = branchTransformer(:,5);
SCUCdata.branchTransformer.K = branchTransformer(:,8);
SCUCdata.branchTransformer.P = branchTransformer(:,7)/standardP;


%%%  24ʱ�θ������ݼ���ת����
% fgetl(filedata);     %����ע��
fgetl(filedata);     %����ע��
checkZero = fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    load(k,:) = fscanf(filedata,'%f',[1,7]);
    k         = k+1;
    checkZero = fscanf(filedata,'%d',1);
end
SCUCdata.totalLoad.PD_T = load(:,1)/standardP;    %%%  ��ʾ T ʱ�θ���ʱ�ε��ܸ���
SCUCdata.totalLoad.QD_T = load(:,2)/standardP;
R_T  = load(:,5:7)/standardP;
SCUCdata.totalLoad.R_T  = sum(R_T,2);                  %% ���ڸ��ɵĲ�������̫�˽⡣������????   ���Ϊ����������Ϊ���������ת����Լ��
SCUCdata.totalLoad.T = 24;

%%%  �ڵ㸺������
% fgetl(filedata);     %����ע��
fgetl(filedata);     %����ע��
checkZero=fscanf(filedata,'%d',1);
k=1;
while checkZero~=0
    nodeSequence(k) = checkZero;
    nodePQ(k,:)     = fscanf(filedata,'%f',[1,2]);
    k               = k+1;
    checkZero       = fscanf(filedata,'%d',1);
end
SCUCdata.busLoad.bus_PDQR = nodeSequence';
SCUCdata.busLoad.node_P   = nodePQ(:,1)/standardP;
SCUCdata.busLoad.node_Q   = nodePQ(:,2)/standardP;

SCUCdata.busLoad.bus_P_factor = SCUCdata.busLoad.node_P/sum(SCUCdata.busLoad.node_P);  %%%  �ڵ㸺������   �����ڵ�ռ���ɵı�������������������ܸ����з��为�ɡ�
SCUCdata.busLoad.bus_Q_factor = SCUCdata.busLoad.node_Q/sum(SCUCdata.busLoad.node_Q);

SCUCdata.busLoad.node_P = SCUCdata.totalLoad.PD_T * SCUCdata.busLoad.bus_P_factor';
SCUCdata.busLoad.node_Q = SCUCdata.totalLoad.QD_T * SCUCdata.busLoad.bus_Q_factor';


N=SCUCdata.units.N;
FileName=strcat('CET\CET_',File);%,'_',num2str(N)
fid_CET=fopen(FileName);
CET= fscanf(fid_CET,'%f',[3,N]);
SCUCdata.a=CET(3,1:N)';
SCUCdata.b=CET(2,1:N)';
SCUCdata.c=CET(1,1:N)';
CET= fscanf(fid_CET,'%f',[1,5]);
SCUCdata.price_buy=CET(1);
SCUCdata.price_sell=CET(2);
SCUCdata.emmission_buy_max=CET(3);
SCUCdata.emmission_sell_max=CET(4);
SCUCdata.E0=CET(5);



end

