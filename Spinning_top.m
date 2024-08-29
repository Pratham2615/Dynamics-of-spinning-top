
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%defining parameters for top
m=1; r=0.05; a=0.005; L=0.1; v=0.01; g=9.81;
Para=[m r a L v g];
tiledlayout(3,2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for finding Θ with different initial values of Θ
theta0=[2 5 10 20 50];
for i=1:length(theta0)
%defining time range
t=[0 0.8];
%defining initial conditions
Y0=[0 theta0(i)*pi/180 0 0 0 50];
%solving the governing equations
[T,Y]=ode45(@(t,Y) top(t,Y,Para),t,Y0);
%extracting Θ from solution of governing equations
theta=Y(:,2)*180/pi;
thetamax=90-atan(r/L)*180/pi;
for j=1:size(theta)
if theta(j)>thetamax
theta(j)=thetamax;
end
end
idx=find(theta==thetamax);
tf(i)=T(idx(1));
%plotting Θ with time
figure;
plot(T,theta)
end
xlabel("Time (s)")
ylabel("Θ (deg)")
title("Variation of Θ with time [ω3(t=0)=50rad/s]")
legend("Θ(0)=2deg","Θ(0)=5deg","Θ(0)=10deg","Θ(0)=20deg","Θ(0)=50deg")
figure;
plot(theta0,tf,"-o")
xlabel("Initial value of Precision Angle (deg)")
ylabel("Time when top touches ground (s)")
title("Variation of touching time with initial precision angle")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for finding Θ with different initial values of ω3
omega30=[50 100 200 500 1000];
for i=1:length(omega30)
%defining time range
t=[0 0.8];
%defining initial conditions
Y0=[0 2*pi/180 0 0 0 omega30(i)];
%solving the governing equations
[T,Y]=ode45(@(t,Y) top(t,Y,Para),t,Y0);
%extracting Θ from solution of governing equations
theta=Y(:,2)*180/pi;
thetamax=90-atan(r/L)*180/pi;
for j=1:size(theta)
if theta(j)>thetamax
theta(j)=thetamax;
end
end
idx=find(theta==thetamax);
tf(i)=T(idx(1));
%plotting Θ with time
figure;
plot(T,theta)
end
xlabel("Time (s)")
ylabel("Θ (deg)")
title("Variation of Θ with time [Θ(t=0)=2deg]")
legend("ω3(0)=50rad/s","ω3(0)=100rad/s","ω3(0)=200rad/s","ω3(0)=500rad/s","ω3(0)=1000rad/s")
figure;
plot(omega30,tf,"-o")
xlabel("Initial value of ω3 (rad/s)")
ylabel("Time when top touches ground (s)")
title("Variation of touching time with initial ω3")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for finding Θ with different values of L/r
Lr=[2 4 6 8 10];
for i=1:length(Lr)
%defining time range
t=[0 1.2];
%defining initial conditions
Para(4)=Lr(i)*Para(2);
Y0=[0 2*pi/180 0 0 0 50];
%solving the governing equations
[T,Y]=ode45(@(t,Y) top(t,Y,Para),t,Y0);
%extracting Θ from solution of governing equations
theta=Y(:,2)*180/pi;
thetamax=90-atan(r/L)*180/pi;
for j=1:size(theta)
if theta(j)>thetamax
theta(j)=thetamax;
end
end
idx=find(theta==thetamax);
%plotting Θ with time
figure;
plot(T,theta)
end
xlabel("Time (s)")
ylabel("Θ (deg)")
title("Variation of Θ with time [Θ(t=0)=2deg and ω3(t=0)=50rad/s]")
legend("L/r=2","L/r=4","L/r=6","L/r=8","L/r=10")
figure;
plot(Lr,tf,"-o")
xlabel("L/r")
ylabel("Time when top touches ground (s)")
title("Variation of touching time with L/r")

function dYdt=top(t,Y,Para)
%defining variables
phi=Y(1); theta=Y(2); psi=Y(3);
omega1=Y(4); omega2=Y(5); omega3=Y(6);
%defining top properties
m=Para(1); r=Para(2); a=Para(3); L=Para(4); v=Para(5); g=Para(6);
%defining moment of inertia about point of contact
I=m*(L^2+0.2*(r^2+a^2));
I0=0.4*m*r^2;
%defining governing equations
if theta>pi/2-atan(r/L)
dthetadt=0;
else
dthetadt=omega1*cos(psi)-omega2*sin(psi);
end
dphidt=(omega1*sin(psi)+omega2*cos(psi))/sin(theta);
dpsidt=omega3-cot(theta)*(omega1*sin(psi)+omega2*cos(psi));
domega1dt=(-v*omega1+m*g*L*cos(psi)*sin(theta)+(I-I0)*omega2*omega3)/I;
domega2dt=(-v*omega2-m*g*L*sin(psi)*sin(theta)-(I-I0)*omega1*omega3)/I;
domega3dt=-v*omega3/I0;
dYdt=[dphidt;dthetadt;dpsidt;domega1dt;domega2dt;domega3dt];
end