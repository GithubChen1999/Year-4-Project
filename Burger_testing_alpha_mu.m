% we are aiming to model the partial differential equation u_t=u_xx=,
% where x is 
%%
clearvars
h=0.01; %mesh size
dt=0.0005; %time mesh
t_0=0;t_max=0.5; %max time interval 
len=2;   %lenth of rod
mu_0=2;alpha=1; %%to change
 
 w=1;epsilon=1/2;
Nr_node=len/h +1; %total numbers of nodes
Node_array=[0:h:len];
N=len/h-1;f=zeros(N,1);
u=zeros(N+2,int32((t_max-t_0)/dt+1));%temp of nodes with different dt
M_t=linspace(t_0,t_max,int32(t_max-t_0)/dt+1);

B=diag((2/3)*h*ones(1,N+2))+ diag(h/6*ones(1,N+1),1) + diag(h/6*ones(1,N+1),-1);
B(1,1)=h/3;B(end,end)=h/3;

C=diag(0.5*ones(1,N+1),+1)+diag(-0.5*ones(1,N+1),-1);
C(1,1)=0.5;C(end,end)=0.5;
C=w*C;
method="av";
%%
% for i = 1:N
%     if i*h>0.2 && i*h<0.5  %initial condition for 凸 shape
%         u(i)=1;
%     end
% end

%  for i =1:len/h+1  % initial value
%      u(i,1)=2*sin((i-1)*h);
%  end

%  for i = 1:1/h   % initial value
%      u(i,1)=2*sin(pi*(i-1)*h);
%  end
    for i = 1:1/h+1   % initial value
        u(i,1)=sin(pi*(i-1)*h);
    end
%%
 %now I use BE to evaluate the time component
time_approximation="RK2";
A_unit=[1 -1;-1,1]./h;

if time_approximation=="RK2"
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
        
        

        C_1=zeros(N+2,1); C_2=zeros(N+2,1);      
        u_1=u(:,i);
        C_1(1)=-2*u_1(1)^2+u_1(1)*u_1(2)+u_1(2)^2;
        C_1(end)=-u_1(end-1)^2-u_1(end)*u_1(end-1)+2*u_1(end)^2;
        for j=2:N+1
            C_1(j)=-2*u_1(j-1)^2-u_1(j-1)*u_1(j)+u_1(j)*u_1(j+1)+2*u_1(j+1)^2;       
        end
        C_1=1/6*C_1;
        k1=-dt*pinv(B)*(C_1+A*u(:,i));
        
        u_2=u(:,i)+k1/2;
        C_2(1)=-2*u_2(1)^2+u_2(1)*u_2(2)+u_2(2)^2;
        C_2(end)=-u_2(end-1)^2-u_2(end)*u_2(end-1)+2*u_2(end)^2;
        for j=2:N+1
            C_2(j)=-2*u_2(j-1)^2-u_2(j-1)*u_2(j)+u_2(j)*u_2(j+1)+2*u_2(j+1)^2;       
        end
        C_2=1/6*C_2;  
        k2=-dt*pinv(B)*(C_2+A*(u(:,i)+k1/2));
        u(:,i+1)=u(:,i)+k2;
      end
end
%%
%  save('saveE.mat',"uE");
% figure 
% plot(Node_array,u(:,0.25*t_max/dt+1))
%  ylim([-1 2])
%  xlim([0 len])
%  title("mesh =" + h  + " dt =" + dt,"viscosity =" +method + " time =" + 0.25*t_max/+1) %+" m="+mu_0+"  alpha="+alpha
%  
% figure 
% plot(Node_array,u(:,0.5*t_max/dt+1))
%  ylim([-1 2])
%  xlim([0 len])
%  title("mesh =" + h  + " dt =" + dt,"viscosity =" +method + " time =" + 0.5*t_max/+1) %+" m="+mu_0+"  alpha="+alpha
%  
%  figure 
% plot(Node_array,u(:,0.75*t_max/dt+1))
%  ylim([-1 2])
%  xlim([0 len])
%  title("mesh =" + h  + " dt =" + dt,"viscosity =" +method + " time =" + 0.75*t_max/+1) %+" m="+mu_0+"  alpha="+alpha
 

 figure 
plot(Node_array,u(:,1*t_max/dt+1))
 ylim([-1 2])
 xlim([0 len])
 title("mesh =" + h  + " dt =" + dt,"viscosity =" +method + " time =" + 1*t_max/+1 +" m="+mu_0+"  alpha="+alpha)
%%
%  figure
%  surf(M_t,Node_array,uG_2)
%  zlim([-1 2])
%  xlabel('time')
%  ylabel('nodes')
%  title("FEM with mesh =" + h  + " timestep =" + dt)
%  hold on
 
 %%

 %figure
 %real value
%  Node=[0:0.01:1];
%  time=linspace(0,0.5,200)';
%  Temp=exp((-pi^2).*time)*sin(pi.*Node);
%  surf(time,Node,Temp')
%  xlabel('time')
%  ylabel('nodes')
%  title("real solution")

 %%
% test = VideoWriter('test.avi');
% open(test)
% figure
%  for i =1 : (t_max-t_0)/dt
%      U_t=u(:,i);
%      plot(Node_array,U_t)
%      ylim([-1 2])
%      xlim([0 len])
%       xlabel('node')
%     ylabel('Temp')
%     title("Burger's equation with mesh =" + h  + " timestep =" + dt,"viscosity =" +method + " time =" + i*dt)
%      F=getframe(gcf);
%      writeVideo(test,F)
%  end
%  close(test)


