% we are aiming to model the partial differential equation u_t=u_xx=,
% where x is 
%%
clearvars
t_0=0;t_max=1; %max time interval 
testh=[ 0.005 0.01 0.02 0.05 0.1];
testdt=[];
errorarray1=[];
errorarray2=[];
len=2; %lenth of rod
method="no";
%%
for v=testh
    h=v; %mesh size
    dt=v/10; %time mesh
    testdt=[testdt,dt];
    mu_0=1/10;epsilon=1/2;alpha=1;
    Nr_node=len/h +1; %total numbers of nodes
    Node_array=[0:h:len];
    N=len/h-1;f=zeros(N,1);
    u=zeros(N+2,int32((t_max-t_0)/dt+1));%temp of nodes with different dt

    B=diag((2/3)*h*ones(1,N+2))+ diag(h/6*ones(1,N+1),1) + diag(h/6*ones(1,N+1),-1);
    B(1,1)=h/3;B(end,end)=h/3;

    C=diag(0.5*ones(1,N+1),+1)+diag(-0.5*ones(1,N+1),-1);
    C(1,1)=0.5;C(end,end)=0.5;
%%
    for i = 1:1/h+1   % initial value
        u(i,1)=sin(pi*(i-1)*h);
    end

%%
 %now I use RK2 to evaluate the time component

    A_unit=[1 -1;-1,1]./h;
 
    for i =1:int32((t_max-t_0)/dt)
        if method=="av"
            d_xU=zeros(Nr_node-1,1);
            for j =1:Nr_node-1
                 d_xU(j)=abs((u(j+1,i)-u(j,i))/h);
            end
            av=h*mu_0*d_xU;
            MINI=min(1,av).^alpha;
            eps=epsilon*h*MINI;
            A=zeros(N+2);
            for j=1:N+1
                A([j,j+1],[j,j+1])=A([j,j+1],[j,j+1])+eps(j)*A_unit;
            end
        end
        if method=="upwind"
            eps=epsilon*h;
            A=zeros(N+2);
            for j=1:N+1
                A([j,j+1],[j,j+1])=A([j,j+1],[j,j+1])+eps*A_unit;
            end
        end
        if method=="no"
            A=0;
        end
        M=A+C;
        D=dt*M+B;
        u(:,i+1)=pinv(D)*(B*u(:,i));
    end
 
 %%

    %real value
    Node=[0:h:len];
    time=[0:dt:t_max]';
    u1=zeros(N+2,1);
    for j=1/h+1:N+2
      u1(j,1)=sin(pi*(j-1)*h-pi);
    end
%     plot(Node,u(:,end))
%     xlabel("Node")
%     ylim([-0.2 1.2])
%     legend("\Delta t="+v)
%     hold on
    u1_=u1([1/h+2:N+1]);
    u_=u([1/h+2:N+1],end);
    Nr_node_= 1/h-1;    
    diffu=u1_-u_(:,end);
    diff1=sum(abs(diffu))/Nr_node_;
    errorarray1=[errorarray1,diff1];
    diff2=sqrt(sum((diffu).^2))/Nr_node_;
    errorarray2=[errorarray2,diff2];
    
end
 %%

figure
p=polyfit(log(testdt),log(errorarray1),1);
loglog(testdt,errorarray1)
xlabel("log(dt)")
ylabel("log(error)")
title("BE with "+method+" viscosity "+"  Error=1-norm" ," Gradient is "+ p(1))

figure
p=polyfit(log(testdt),log(errorarray2),1);
loglog(testdt,errorarray2)
xlabel("log(dt)")
ylabel("log(error)")
title("BE with "+method+" viscosity "+"  Error=2-norm" ," Gradient is "+ p(1))


