r=51; % Number of rows
c=4*r; % Number of columns
t_final=100;

%Given
q_dot=100*100^2; %[W/m^2]
k=120; % [W/m.C]
rho=8500; % [kg/m^3]
Cp=400; % [J/kg.K]
T_inf=10; % [C]

dx=(4/100)/(c-1); %[m]
dy=(1/100)/(r-1); %[m]
alpha=k/rho/Cp; %[m^2/K]
dt=1/(2*alpha*(1/dx^2+1/dy^2)); %[sec] Stability Criterion
beta=alpha*dt/dx^2;
gamma=alpha*dt/dy^2;

%Exact Solution
Te=zeros(r,1);
Te(:,1)=T_inf;
Tge=zeros(r,1);
Tge(:,1)=Te(:,1);
for y=0:dy:0.01
    ii=2:1:r-1;
    Tge(1,1)=q_dot*y/k+Te(1,1);
    Tge(ii,1)=Tge(1,1)-q_dot*(dy*(ii-1))/k;
end
ExactSolution_MaxTemp=Tge(1)

%Numerical Solution
for i=1:1:r
    j=1:1:c;
    T=zeros(r,c);
    T(1,:)=T_inf;
    T(i,:)=T_inf;
    T(:,1)=T_inf;
    T(:,j)=T_inf;
end
Tg=zeros(r,c);
Tg(1,:)=T(1,:);
Tg(i,:)=T(i,:);
Tg(:,1)=T(:,1);
Tg(:,j)=T(:,j);

for t=0:dt:t_final
    i=2:1:r-1;
    j=2:1:c-1;
    Tg(1,j)=(1-2*beta-2*gamma)*T(1,j)+beta*(T(1,j-1)+T(1,j+1))+2*gamma*T(2,j)+2*q_dot*dt/(rho*Cp*dy); %Top boundary equation
    Tg(i,j)=(1-2*beta-2*gamma)*T(i,j)+beta*(T(i,j-1)+T(i,j+1))+gamma*(T(i+1,j)+T(i-1,j)); %Internal nodes
    if (Tg-T)<1e-5
        break
    end
    T=Tg;
    tf=0:dt:t;
    iter=length(tf);
    Tgmax(iter)=max(Tg(1,:));
end

MaxTemp=Tgmax(end)
FinalTime=t
percentoferror=(Tge(1)-Tgmax(end))/Tge(1)*100

figure;
x=0:dx:0.01;
plot(x,T(:,c/2),x,Tge)
legend('Numerical Method','Exact Solution')
xlabel('Distance from top wall to bottom wall')
ylabel('Temperature Distribution')
figure;
plot(tf,Tgmax)
xlabel('Time(sec)')
title('Temperature distribution along a line running through the center')
ylabel('Maximum Temperture')
title(['t = ' num2str(t)])
figure;
image(T)
title(['t = ' num2str(t)])
colorbar
