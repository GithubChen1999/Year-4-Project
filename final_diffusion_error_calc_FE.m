% we are aiming to model the partial differential equation u_t=u_xx=,
% where x is 
%%
clearvars

t_0=0;t_max=0.2; %max time interval 
len=1;   %lenth of rod
u_0=0; u_end=0;  %Temp at the start and end
w=0;

errorarray=[]
%%
testh=[0.005 0.01 0.02 0.05 0.1];
testdt=1/6*testh.^2;
for v=testh
    h=v; %mesh size
    Nr_node=len/h +1; %total numbers of nodes
    Node_array=[0:h:len];
    N=len/h-1;f=zeros(N,1);
    A =diag(2*ones(1,N)) + diag(-1*ones(1,N-1),1) + diag(-1*ones(1,N-1),-1);
    A=(1/h)*A;
    B=diag((2/3)*h*ones(1,N))+ diag(h/6*ones(1,N-1),1) + diag(h/6*ones(1,N-1),-1);

    dt=1/6*v^2; %time mesh
    u=zeros(N,int32((t_max-t_0)/dt+1));%temp of nodes with different dt
     for i = 1:N   % initial value
         u(i,1)=sin(pi*i*h);
     end
    M_t=linspace(t_0,t_max,int32((t_max-t_0)/dt+1));



    % %now I use FE to evaluate the time component
      for i =1:int32((t_max-t_0)/dt)
         u(:,i+1)=inv(B)*((B-dt*A)*u(:,i));
     end

          
     u=[u_0*ones(1,int32((t_max-t_0)/dt+1));u];
     u=[u;u_end*ones(1,int32((t_max-t_0)/dt+1))];
%%
    Node=[0:h:1];
    time=[0:dt:t_max]';
    Temp=exp((-pi^2).*time)*sin(pi.*Node);
    u1=Temp';t=time';n=Node';

     diff=sqrt(sum(abs(u1(:,end)-u(:,end)).^2));

    errorarray=[errorarray,diff/Nr_node];
     
     
end
%%
figure
p=polyfit(log(testdt),log(errorarray),1);
loglog(testdt,errorarray)
xlabel("log(dt)")
ylabel("log(error)")
title("FE: Gradient is "+ p(1) )