function [x,fobj,total_time]=function_PSO(model,dataUC,timelimit)
% Define the details of the binary optimization problem
tic
nVar = length(model.ctype);
ub =model.ub;
lb = model.lb;
IntCon=find(model.ctype=='B');
num_int=length(IntCon);

% Define the PSO's paramters
noP = 50;
maxItert = 10000;
wMax = 0.9;
wMin = 0.4;
c1_min=1.5;
c1_max=2.2
c2_min=1.5;
c2_max=2.2;
count_Best_history=0;

% fai=c1+c2;
% w=2/abs(2-fai-(fai^2-4*fai)^0.5);
vMax = repmat((ub - lb) .* 0.2,1,noP);
vMin  = -vMax;
r_eq=1e-4;
total_time=0;

%SR and min_on_off coefficient
Alpha = dataUC.alpha;                           %火电机组发电函数系数Alpha--N*1矩阵
Beta = dataUC.beta;                             %火电机组发电函数系数Beta--N*1矩阵
Gama = dataUC.gamma;                            %火电机组发电函数系数Gama--N*1矩阵
ThPimin = dataUC.p_low;                         %火电机组发电功率下界--N*1矩阵
ThPimax = dataUC.p_up;                          %火电机组发电功率上界--N*1矩阵
Alpha_wan = Alpha + Beta .* ThPimin + Gama .* ThPimin .* ThPimin;           %Alpha弯弯--N*1矩阵
Beta_wan = (ThPimax - ThPimin) .* (Beta + 2 * Gama .* ThPimin);             %Beta弯弯--N*1矩阵
Gama_wan = Gama .* (ThPimax - ThPimin) .* (ThPimax - ThPimin);

SR=dataUC.spin;
Ton=dataUC.time_min_on;
Toff=dataUC.time_min_off;
Tint=dataUC.time_on_off_ini;
T=dataUC.T;
N=dataUC.N;
a=Alpha_wan;
b=Beta_wan;
c=Gama_wan;
Pmax_i=dataUC.p_up;
u0 = full(spones(dataUC.p_initial)); 
Ui = max(0,min(ones(N,1) * T,u0 .* (Ton - Toff)));                    %--N*1矩阵
Li = max(0,min(ones(N,1) * T,(ones(N,1) - u0) .*  (Toff + Tint)));     %--N*1矩阵


min_per_unit_cost=(2*(c.*a).^0.5+b)./Pmax_i;
Plist=[];
while(any(min_per_unit_cost))
    max_P=max(min_per_unit_cost);
    index_max_P=find(min_per_unit_cost==max_P);
    Plist=[index_max_P;Plist];
    min_per_unit_cost(index_max_P)=0;
end

segma=sparse(T,1);
uit=sparse(N,T);
for i=1:N
    if(Ui(i) + Li(i) >= 1)
        for t=1:Ui(i)+Li(i)
            uit(i,t)=1;
        end
    end
end
for t=1:T
    for i=1:N
        while segma(t)<SR(t)
            uit(Plist(i),t)=1;
            segma(t)=segma(t)+Pmax_i(Plist(i));
        end
    end
end

%step 1
t_on=sparse(N,T);
t_off=sparse(N,T);
sit=sparse(N,T);
for i=1:N
    for t=1:T
        if(t==1) %t=1
            if(uit(i,t)==1)
                if(u0(i)==1)
                    t_on(i,t)=Tint(i)+1;
                else
                    t_on(i,t)=1;
                end
            else
                if(u0(i)==0)
                    t_off(i,t)=Tint(i)-1;
                else
                    t_off(i,t)=-1;
                end
            end
        else %t>1
            if(uit(i,t)==1)
                t_on(i,t)=t_on(i,t-1)+1;
            else
                t_off(i,t)=t_off(i,t-1)-1;
            end
        end
    end
end


for t=2:T
    for i=1:N
        if (uit(i,t)==0&uit(i,t-1)==1&t_on(i,t-1)>=Ton(i)-t_on(i,t-1))
            uit(i,t)=1;
        else
            uit(i,t-1)=0;
        end
        if(uit(i,t-1)==0&uit(i,t)==1&abs(t_off(t-1))<Toff)
            q=abs(t_off(t-1));
            if(t-q>0)
                uit(i,t-q:t-1)=1;
            else
                uit(i,1:t-1)=1;
            end
        end
        if(uit(i,t)==1)
            t_on(i,t)=t_on(i,t-1)+1;
            t_off(i,t)=0;
        else
            t_off(i,t)=t_off(i,t-1)-1;
            t_on(i,t)=0;
        end
    end
end

for i=1:N
    for t=2:T
        if(uit(i,t)==1&uit(i,t-1)==0)
            sit(i,t)=1;
        end
    end
end



var_int=[reshape(uit',[N*T,1]);reshape(sit',[N*T,1])];

index_lb_in_inf=find(lb~=-inf);
index_ub_in_inf=find(ub~=inf);
constraint_lb=sparse(1:length(index_lb_in_inf),index_lb_in_inf,-1,length(index_lb_in_inf),nVar);
constraint_ub=sparse(1:length(index_ub_in_inf),index_ub_in_inf,1,length(index_ub_in_inf),nVar);

A=[model.Aineq;model.Aeq;-1*model.Aeq;constraint_lb;constraint_ub];
b=[model.bineq;model.beq+r_eq*ones(length(model.beq),1);r_eq*ones(length(model.beq),1)-1*model.beq;-1*lb(index_lb_in_inf);ub(index_ub_in_inf)];


% Initialize the particles
index_lb_inf=find(lb==-inf);
index_ub_inf=find(ub==inf);
lb(index_lb_inf)=-noP*nVar;
ub(index_ub_inf)=noP*nVar;

X_Swarm_Particles=repmat((ub-lb),1,noP).*rand(nVar,noP)+repmat(lb,1,noP);
X_Swarm_Particles(IntCon,:)=repmat(var_int,1,noP);
V_Swarm_Particles=sparse(nVar,noP);
PBEST_X_Swarm_Particles=sparse(nVar,noP);
PBEST_Obj_Swarm_Particles=inf*ones(noP,1);
GBEST_X_Swarm_Particles=sparse(nVar,1);
GBEST_Obj_Swarm_Particles=inf;

time=toc;
total_time=total_time+time;

% Main loop
for iter = 1 : maxItert
    tic
    % Calcualte the objective value
    currentX = X_Swarm_Particles;
    g=A*currentX-repmat(b,1,noP);
    for i=1:noP
        g_noncon(i)=currentX(:,i)'*model.Q*currentX(:,i)+currentX(:,i)'*model.l;
    end
    g=[g;g_noncon];
    unfeasiblity_con=find(g>r_eq);
    q=max(sparse(size(g,1),noP),g);
    [index_seita1_row,index_seita1_col]=find(q<0.1);
    [index_seita2_row,index_seita2_col]=find(0.1<=q&q<=1);
    [index_seita3_row,index_seita3_col]=find(q>1);
    seita1=sparse(index_seita1_row,index_seita1_col,10,size(g,1),noP);
    seita2=sparse(index_seita2_row,index_seita2_col,20,size(g,1),noP);
    seita3=sparse(index_seita3_row,index_seita3_col,100,size(g,1),noP);
    seita=seita1+seita2+seita3;
    [index_gamma1_row,index_gamma1_col]=find(q<1);
    [index_gamma2_row,index_gamma2_col]=find(q>=1);
    gamma1=sparse(index_gamma1_row,index_gamma1_col,1,size(g,1),noP);
    gamma2=sparse(index_gamma2_row,index_gamma2_col,2,size(g,1),noP);
    gamma=gamma1+gamma2;
    C=iter^(0.5);
    Obj_f=model.f'*currentX;
    Obj_Swarm_Particles = (model.f'*currentX+C*sum(seita.*q.*exp(gamma)))';
    
    index_PBEST_X=find(Obj_Swarm_Particles-PBEST_Obj_Swarm_Particles<0&Obj_f'>0);
    PBEST_Obj_Swarm_Particles(index_PBEST_X)=Obj_Swarm_Particles(index_PBEST_X);
    PBEST_X_Swarm_Particles(:,index_PBEST_X)=X_Swarm_Particles(:,index_PBEST_X);
    
    index_PWORST_X=find(Obj_Swarm_Particles-PBEST_Obj_Swarm_Particles>0);
    if(isempty(index_PWORST_X)==0)
    PWORST_Obj_Swarm_Particles(index_PWORST_X)=Obj_Swarm_Particles(index_PWORST_X);
    PWORST_X_Swarm_Particles(:,index_PWORST_X)=X_Swarm_Particles(:,index_PWORST_X);
    else
        PWORST_X_Swarm_Particles=X_Swarm_Particles;
    end
    index_GBEST_X=find(Obj_Swarm_Particles<GBEST_Obj_Swarm_Particles);
    if(isempty(index_GBEST_X)==0)
        GBEST_Obj_Swarm_Particles=min(Obj_Swarm_Particles(index_GBEST_X));
        index_GBEST=find(Obj_Swarm_Particles==GBEST_Obj_Swarm_Particles);
        GBEST_X_Swarm_Particles=X_Swarm_Particles(:,index_GBEST);
    end
    
    % Update the X and V vectors
         w = wMin+((maxItert-iter) * ((wMax - wMin)) / maxItert);
         w_new=wMin+w*rand();
         c1_b=c1_min+(c1_max-(c1_max-c1_min)*(iter/maxItert))*rand();
         c1_g=c1_min+(c1_max-(c1_max-c1_min)*(iter/maxItert))*rand();
         c2=c2_min+(c2_max-(c2_max-c2_min)*(iter/maxItert))*rand();
    V_Swarm_Particles=w_new*V_Swarm_Particles...
        +c1_b*rand(nVar,noP).*(PBEST_X_Swarm_Particles-X_Swarm_Particles)...
        +c1_g*rand(nVar,noP).*(X_Swarm_Particles-PWORST_X_Swarm_Particles)...
        +c2*rand(nVar,noP).*(repmat(GBEST_X_Swarm_Particles,1,noP)-X_Swarm_Particles);
    miu=wMin-exp((-w_new)/wMax);
    if(miu>rand()||count_Best_history>=5)
        V_Swarm_Particles=vMax.*rand(nVar,noP);
        count_Best_history=0;
    end
    [index_Max_row,index1_Max_col] = find(V_Swarm_Particles > vMax);             % Check velocities
    [index_Min_row,index1_Min_col] = find(V_Swarm_Particles < vMin);
    for i=1:length(index_Max_row)
        V_Swarm_Particles(index_Max_row(i),index1_Max_col(i)) = vMax(index_Max_row(i),index1_Max_col(i));
    end
    for i=1:length(index_Min_row)
        V_Swarm_Particles(index_Min_row(i),index1_Min_col(i)) = vMin(index_Min_row(i),index1_Min_col(i));
    end
    X_Swarm_Particles=X_Swarm_Particles+V_Swarm_Particles;
    for i=1:noP
        index_x_lb=find(X_Swarm_Particles(:,i)-lb<0);
        index_x_ub=find(X_Swarm_Particles(:,i)-ub>0);
        X_Swarm_Particles(index_x_lb,i)=(ub(index_x_lb)-lb(index_x_lb))*rand()+lb(index_x_lb);
        X_Swarm_Particles(index_x_ub,i)=(ub(index_x_ub)-lb(index_x_ub))*rand()+lb(index_x_ub);
    end
    s = 1 ./ (1 + exp(V_Swarm_Particles));
    
    r = rand(num_int,noP);
    [index_int_row,index_int_col]=find(r < s(IntCon,:));
    for i=1:length(index_int_row)
        X_Swarm_Particles(index_int_row(i),index_int_col(i))=0;
    end
    
    [index_int_row,index_int_col]=find(r >= s(IntCon,:));
    for i=1:length(index_int_row)
        X_Swarm_Particles(index_int_row(i),index_int_col(i))=1;
    end


    for k=1:noP
        t_on=sparse(N,T);
        t_off=sparse(N,T);
        sit=sparse(N,T);
        uit=reshape(X_Swarm_Particles(1:N*T,k),[T,N])';
        for i=1:N
            if(Ui(i) + Li(i) >= 1)
                for t=1:Ui(i)+Li(i)
                    uit(i,t)=1;
                end
            end
        end
    for i=1:N
        for t=1:T
            if(t==1) %t=1
                if(uit(i,t)==1)
                    if(u0(i)==1)
                        t_on(i,t)=Tint(i)+1;
                    else
                        t_on(i,t)=1;
                    end
                else
                    if(u0(i)==0)
                        t_off(i,t)=Tint(i)-1;
                    else
                        t_off(i,t)=-1;
                    end
                end
            else %t>1
                if(uit(i,t)==1)
                    t_on(i,t)=t_on(i,t-1)+1;
                else
                    t_off(i,t)=t_off(i,t-1)-1;
                end
            end
        end
    end
    
    for t=2:T
        for i=1:N
            if (uit(i,t)==0&uit(i,t-1)==1&t_on(i,t-1)>=Ton(i)-t_on(i,t-1))
                uit(i,t)=1;
            else
                uit(i,t-1)=0;
            end
            if(uit(i,t-1)==0&uit(i,t)==1&abs(t_off(t-1))<Toff)
                q=abs(t_off(t-1));
                if(t-q>0)
                    uit(i,t-q:t-1)=1;
                else
                    uit(i,1:t-1)=1;
                end
            end
            if(uit(i,t)==1)
                t_on(i,t)=t_on(i,t-1)+1;
                t_off(i,t)=0;
            else
                t_off(i,t)=t_off(i,t-1)-1;
                t_on(i,t)=0;
            end
        end
    end
    
    for i=1:N
        for t=2:T
            if(uit(i,t)==1&uit(i,t-1)==0)
                sit(i,t)=1;
            end
        end
    end
    var_int=[reshape(uit',[N*T,1]);reshape(sit',[N*T,1])];
    X_Swarm_Particles(IntCon,k)=var_int;
    end
    
    if(iter>1&model.f'*GBEST_X_Swarm_Particles==Obj_Best_history)
        count_Best_history=count_Best_history+1;
    end
    if(iter>1&model.f'*GBEST_X_Swarm_Particles~=Obj_Best_history)
        count_Best_history=0;
    end
    time=toc;
    total_time=total_time+time;
    outmsg = ['Iteration# ', num2str(iter) , ' Obj = ' , num2str(model.f'*GBEST_X_Swarm_Particles),' Time=',num2str(total_time),'  unfeasibility  ',num2str(length(unfeasiblity_con))];
    disp(outmsg);
    Obj_Best_history=model.f'*GBEST_X_Swarm_Particles;
    if(total_time>timelimit)
        x=GBEST_X_Swarm_Particles;
        fobj=model.f'*GBEST_X_Swarm_Particles;
        
        return;
    end
end

x=GBEST_X_Swarm_Particles;
fobj=model.f'*GBEST_X_Swarm_Particles;
end


