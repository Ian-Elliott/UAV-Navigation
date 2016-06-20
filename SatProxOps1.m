function SatProxOps1
%%Mixed integer linear programming of trajectory in 3 dimensional space
% with obstacle avoidance. Bounded by initial and final postion and 
% velocity, max thrust, and total time.
% 
% Accepts values for total time (tf), time step (dt), inital and final
% states (x0,v0,xf,vf), max thrust for each axis (uxmax), and mass (m).
% Also accepts the locations of multiple objects to avoid collsion.
clc,close all

%% Inputs
% Simulation time
tf = 100; %s
dt = 1; %s

% Orbital parameters
a = 10e6; %m
n = a2n(a); %rad/s

% Initial state
x0 = -50; %m
y0 = 3;
z0 = 0;
vx0 = 0; %m/s
vy0 = 0;
vz0 = 0;

% Final state
xf = 50; %m
yf = 0;
zf = 0;
vxf = 0; %m/s
vyf = 0;
vzf = 0;

% Vehicle parameters
uxmax = 1; %N
uymax = uxmax;
uzmax = 1;
m = 12; %kg

% Obstacles bounds
% xbmin = [-5 -2 -5.5 5]; %m
% ybmin = [-0.5 -2 -5 -5];
% zbmin = [-0.5 -2 -5 -5];
% xbmax = [5 2 -5 5.5];
% ybmax = [0.5 2 5 5];
% zbmax = [0.5 2 5 5];
xbmin = [-5,-0.5,-0.5]; %m
ybmin = [-5,5,-30];
zbmin = [-5,-2.5,-2.5];
xbmax = [5,0.5,0.5];
ybmax = [5,30,-5];
zbmax = [5,2.5,2.5];

% Safety buffer
d = 0.5; %m

%% Pre-optimization setup
t = 0:dt:tf;
Nsim = length(t)-1;
Nvar = 6*Nsim;
if exist('xbmin','var')
    NObj = length(xbmin);
else
    NObj = 0;
end
Nbi = NObj*Nvar;

% Time constants
alpha = dt/m;
beta = dt^2/m;

% Parameter bounds, function coefficients, integer contraints
lb = zeros(Nvar+Nbi,1);
ub = ones(Nvar+Nbi,1);
f = zeros(Nvar+Nbi,1);
for i = 1:6:Nvar
    ub(i:i+1) = [uxmax;uxmax];
    ub(i+2:i+3) = [uymax;uymax];
    ub(i+4:i+5) = [uzmax;uzmax];
    f(i:i+5) = dt*ones(6,1);
end
intcon = (Nvar+1:Nvar+Nbi);

% Acceleration from HCW
[ax,ay,az] = HCW([x0,y0,z0],[vx0,vy0,vz0],n);
% ax = 0;ay = 0;az = 0;

% Boundary conditions equality constraints
Aeq = zeros(6,Nvar+Nbi);
beq = zeros(6,1);
j = 1;
for i = 1:Nsim
    Aeq(1,j) = beta*(Nsim-i);
    Aeq(1,j+1) = -beta*(Nsim-i);
    Aeq(2,j+2) = beta*(Nsim-i);
    Aeq(2,j+3) = -beta*(Nsim-i);
    Aeq(3,j+4) = beta*(Nsim-i);
    Aeq(3,j+5) = -beta*(Nsim-i);
    Aeq(4,j) = alpha;
    Aeq(4,j+1) = -alpha;
    Aeq(5,j+2) = alpha;
    Aeq(5,j+3) = -alpha;
    Aeq(6,j+4) = alpha;
    Aeq(6,j+5) = -alpha;
        
    j = j+6;
end
c = (Nsim-1):-1:0;
beq(1) = (xf-x0)-Nsim*dt*vx0-dt^2*sum(ax*c);
beq(2) = (yf-y0)-Nsim*dt*vy0-dt^2*sum(ay*c);
beq(3) = (zf-z0)-Nsim*dt*vz0-dt^2*sum(az*c);
beq(4) = (vxf-vx0)-dt*sum(ax);
beq(5) = (vyf-vy0)-dt*sum(ay);
beq(6) = (vzf-vz0)-dt*sum(az);

% Add safety buffer
if NObj > 0
    xbmind = xbmin-d;
    ybmind = ybmin-d;
    zbmind = zbmin-d;
    xbmaxd = xbmax+d;
    ybmaxd = ybmax+d;
    zbmaxd = zbmax+d;
end

% Add obstacle inequality contraints 
A = [];
b = [];
for N = 1:NObj
    [A,b] = NewObstacle(N,A,b,[x0,y0,z0],[vx0,vy0,vz0],n,dt,m,Nsim,Nbi,...
        xbmind(N),xbmaxd(N),ybmind(N),ybmaxd(N),zbmind(N),zbmaxd(N));
end

%% MILP Optimization
sA = size(A);
fprintf('A: %d x %d\n',sA(1),sA(2))
fprintf('u: %d\n',Nvar)
fprintf('int: %d\n',Nbi')
options = optimoptions(@intlinprog,'Display','iter');
[u,fval] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);
fprintf('J: %.3f\n',fval)

%% Post-process
ux = u(1:6:Nvar)-u(2:6:Nvar);
uy = u(3:6:Nvar)-u(4:6:Nvar);
uz = u(5:6:Nvar)-u(6:6:Nvar);

% Binary values
bi = u(intcon);
bimat = zeros(length(bi)/6,6);
i = 1;
for j = 1:6:length(bi)
   bimat(i,:) = round(bi(j:j+5));
   i = i+1;
end

% Coordinates & velocity
x = zeros(Nsim+1,1);y = zeros(Nsim+1,1);z = zeros(Nsim+1,1);
vx = zeros(Nsim+1,1);vy = zeros(Nsim+1,1);vz = zeros(Nsim+1,1);
x(1) = x0;y(1) = y0;z(1) = z0;
vx(1) = vx0;vy(1) = vy0;vz(1) = vz0;
for i = 2:length(ux)
    vx(i) = vx(i-1)+(ux(i-1)/m+ax)*dt;
    x(i) = x(i-1)+vx(i-1)*dt;
    vy(i) = vy(i-1)+(uy(i-1)/m+ay)*dt;
    y(i) = y(i-1)+vy(i-1)*dt;
    vz(i) = vz(i-1)+(uz(i-1)/m+az)*dt;
    z(i) = z(i-1)+vz(i-1)*dt;
end
x(end) = xf;y(end) = yf;z(end) = zf;
vx(end) = vxf;vy(end) = vyf;vz(end) = vzf;
ux = [ux;0];uy = [uy;0];uz = [uz;0];

%% Plots
% Trajectory
figure(1)
hold on
p1 = plot3(x,y,z,'bo','linewidth',2);
p2 = plot3(x0,y0,z0,'k^','linewidth',2,'markersize',10);
p3 = plot3(xf,yf,zf,'kx','linewidth',2,'markersize',10);
p4 = quiver3(x,y,z,-ux,-uy,-uz,0.3,'r','linewidth',2);
hold off
for N = 1:NObj
    PlotObstacle(xbmin(N),xbmax(N),ybmin(N),ybmax(N),zbmin(N),zbmax(N));
end
grid on
xlabel('x, m')
ylabel('y, m')
zlabel('z, m')
title('Fuel Optimal Trajectory')
legend([p1,p2,p3,p4],{'Trajectory','x0','xf','Thrust'},'location','BestOutside')
axis('equal')
camva(9)
view(-45,15)

% Thrust
figure(2)
subplot(3,1,1)
hold on
plot(t,ux/uxmax,'-b','linewidth',2)
plot([0 tf],[1 1],'--k','linewidth',2)
plot([0 tf],[-1 -1],'--k','linewidth',2)
grid on
axis([0 tf -1 1])
title('Control Signals')
ylabel('ux')

subplot(3,1,2)
hold on
plot(t,uy/uymax,'-r','linewidth',2)
plot([0 tf],[1 1],'--k','linewidth',2)
plot([0 tf],[-1 -1],'--k','linewidth',2)
axis([0 tf -1 1])
grid on
ylabel('uy')

subplot(3,1,3)
hold on
plot(t,uz/uzmax,'-g','linewidth',2)
plot([0 tf],[1 1],'--k','linewidth',2)
plot([0 tf],[-1 -1],'--k','linewidth',2)
axis([0 tf -1 1])
grid on
xlabel('Time, s')
ylabel('uz')

% Velocity
figure(3)
subplot(3,1,1)
plot(t,vx,'-b','linewidth',2)
grid on
ylabel('Vx')
title('Velocity')

subplot(3,1,2)
plot(t,vy,'-r','linewidth',2)
grid on
ylabel('Vy')

subplot(3,1,3)
plot(t,vz,'-g','linewidth',2)
grid on
xlabel('Time, s')
ylabel('Vz')
end

function [ax,ay,az] = HCW(p0,v0,n)
% Intial position & velocity
x0 = p0(1);y0 = p0(2);z0 = p0(3);
vx0 = v0(1);vy0 = v0(2);vz0 = v0(3);
% HCW acceleration
ax = 3*x0*n^2+2*vy0*n;
ay = -2*n*vx0;
az = -n^2*z0;
end

function [Anew,bnew] = NewObstacle(N,Aold,bold,p0,v0,n,dt,m,Nsim,Nbi,xbmind,xbmaxd,ybmind,ybmaxd,zbmind,zbmaxd)
% Obstacles inequalities
x0 = p0(1);y0 = p0(2);z0 = p0(3);
vx0 = v0(1);vy0 = v0(2);vz0 = v0(3);
Nvar = 6*Nsim;
beta = dt^2/m;
M = 1e6;

% Acceleration from HCW
[ax,ay,az] = HCW([x0,y0,z0],[vx0,vy0,vz0],n);

% Add inequality contraints
i = 1;
A = zeros(Nsim*7,Nvar+Nbi);
b = zeros(Nsim*7,1);
for k = 1:7:(7*Nsim)
    j = 1;
    for n = 1:i
        % Postive bounds
        A(k,j) = beta*(i-n);
        A(k,j+1) = -beta*(i-n);
        A(k+1,j+2) = beta*(i-n);
        A(k+1,j+3) = -beta*(i-n);
        A(k+2,j+4) = beta*(i-n);
        A(k+2,j+5) = -beta*(i-n);
        
        % Negative bounds
        A(k+3,j) = -beta*(i-n);
        A(k+3,j+1) = beta*(i-n);
        A(k+4,j+2) = -beta*(i-n);
        A(k+4,j+3) = beta*(i-n);
        A(k+5,j+4) = -beta*(i-n);
        A(k+5,j+5) = beta*(i-n);
        j = j+6;
    end
    
    % Binary variables
    j = (Nvar+1)+Nvar*(N-1)+6*(i-1);
    A(k,j) = -M;
    A(k+1,j+1) = -M;
    A(k+2,j+2) = -M;
    A(k+3,j+3) = -M;
    A(k+4,j+4) = -M;
    A(k+5,j+5) = -M;
    A(k+6,j:j+5) = ones(1,6);
    
    c = (i-1):-1:0;
    b(k) = (xbmind-x0)-i*dt*vx0-dt^2*sum(ax*c);
    b(k+1) = (ybmind-y0)-i*dt*vy0-dt^2*sum(ay*c);
    b(k+2) = (zbmind-z0)-i*dt*vz0-dt^2*sum(az*c);
    b(k+3) = -((xbmaxd-x0)-i*dt*vx0-dt^2*sum(ax*c));
    b(k+4) = -((ybmaxd-y0)-i*dt*vy0-dt^2*sum(ay*c));
    b(k+5) = -((zbmaxd-z0)-i*dt*vz0-dt^2*sum(az*c));
    b(k+6) = 5;
    
    i = i+1;
end

% Update matrices
Anew = [Aold;A];
bnew = [bold;b];
end

function PlotObstacle(xbmin,xbmax,ybmin,ybmax,zbmin,zbmax)
% Plot building
xb = [xbmin xbmax xbmax xbmin xbmin];
yb = [ybmin ybmin ybmax ybmax ybmin];
zbl = [zbmin zbmin zbmin zbmin zbmin];
zbu = [zbmax zbmax zbmax zbmax zbmax];

figure(1)
hold on
plot3(xb,yb,zbl,'k','linewidth',2)
plot3(xb,yb,zbu,'k','linewidth',2)
for i = 1:4
    plot3([xb(i) xb(i)],[yb(i) yb(i)],[zbmin,zbmax],'k','linewidth',2)
end
hold off
axis('equal')
camva(9)
view(-45,15)
end

function n = a2n(a)
mu = 3.986004418e14;
n = sqrt(mu/a^3);
end