%1-D
%Given
q_dot=100*100^2; %[W/m^2]Heat Flux on top surface
k=120; % [W/m.C]
rho=8500; % [kg/m^3]
Cp=400; % [J/kg.K] Heat Capacity
T_inf=10; % [C]Temperature maintained at 10 degrees Celsius for left, right and bottom surfaces.

r=51; % Number of rows

dy=(1/100)/(r-1); %[m] Mesh size x-direction
alpha=k/rho/Cp; %[m^2/K]
dt=dy^2/(2*alpha); %[sec] Time step is estimated using Stability Criterion
beta=alpha*dt/dy^2;

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

%Numerical Method
T=zeros(r,1);
T(:,1)=T_inf;
Tg=zeros(r,1);
Tg(:,1)=T(:,1);
t_final=100;
for t=0:dt:t_final
    i=2:1:r-1;
    Tg(1,1)=(1-beta)*T(1,1)+beta*T(2,1)+q_dot*dt/(rho*Cp*dy);
    Tg(i,1)=(1-2*beta)*T(i,1)+beta*(T(i-1,1)+T(i+1,1));
    if (Tg-T)<1e-5
        break
    end
    T=Tg;
    tf=0:dt:t;
    iter=length(tf);
    Tgmax(iter)=max(Tg(1,1));
end

NumericalMaxTemp=Tgmax(end)
FinalTime=t
percentoferror=(Tge(1)-Tgmax(end))/Tge(1)*100

%Comparing Numerical solutions to Exact Solutions
y=0:dy:0.01;
figure;
plot(y,T(:,1),y,Tge(:,1))
legend('Numerical Method','Exact Solution')
xlabel('Distance from top wall to bottom wall')
ylabel('Temperature Distribution')
