%% Coupled Processes CIE4365
% LOTKA - VOLTERRA
%
% Script introducing Numerical Approximation techniques for solving
% Ordinary Differential Equations

%clear memory and close all open figure windows
clear
close all

difftime = eps; % this is assumed to be zero (equality for time)
convcrit = 1e-8; % test value for assessing convergence 

%% Discretization
%Initialize Solver (discretization of space and time)
%No spatial discretization required...

% Time Discretization
deltmax = 0.005;  %max time step;
outtime = (0:.25:250); %times for storing model output
t = 0; %initial time
tend = outtime(end);


%% Initial Model States and model parameters
%Define Parameters for the Lotka Volterra model (see Wikipedia
%and other internet sources

Par.a = 1;
Par.b = 0.2;
Par.c = 0.1;
Par.d = 0.01;

% initial values can be identified from equilibrium state
% equilibrium (y=a/b; x = c/d)
% yini = [Par.c/Par.d Par.a/Par.b]
yini = [11;3];

%% Explicit EULER
%Using explicit Euler

disp('explicit Euler');
tic
%Initialize output matrices
%store initial values in output matricese
T1(1) = 0;

%Note that yini is a column vector, output is stored rows (see ode help)
Y1(1,1:2) = yini';  
y = yini;
nout = 1;

%Run model over time period
while abs(t-tend) > difftime %this approach is best to test if doubles are equal
   %Calc Rates
   dy = LotkaVolterraTHe(t,y,Par);
   
   %Time step should be small enough to prevent negative states from
   %occuring
   %Check delt
   dttest = abs(0.1*y./dy.*(dy<0))+(dy>=0); %we do not allow negative values
   dtout = outtime(nout+1)-t; %we do not want to miss an output time
   delt = min([dttest(:)' deltmax dtout]);
   
   %Update states
   y = y + dy.*delt;
   t = t + delt;
   %Update output matrix
   if abs(t-outtime(nout+1)) < difftime
      nout = nout+1;
      T1(nout) = t;
      Y1(nout,1) = y(1); Y1(nout,2) = y(2);
   end
end
toc

%% Predictor Corrector
% Approach based on Predictor corrector from Wikipedia
disp ('Euler Predictor Corrector');
deltmax = 0.1;
tic
y = yini; %initial states;
t = 0;
%store initial values in output matricese
T2(1) = 0;
Y2(1,1:2) = yini'; 
nout = 1;

while abs(t-tend) > difftime
   %Calc Rates
   dy = LotkaVolterraTHe(t,y,Par);
   
   %Check delt
   dttest = abs(0.1*y./dy.*(dy<0))+(dy>=0);
   dtout = outtime(nout+1)-t;
   delt = min([dttest(:)' deltmax dtout]);
   
   %iteration with predictor corrector
   %start with Euler step
   converged = false;
   yn = y + dy.*delt;
   while ~converged
       dyn = LotkaVolterraTHe(t+delt,yn,Par);
       ynn = y+delt./2.*(dy+dyn);
       if abs(yn-ynn)<convcrit
           converged=true;
       else
           yn=ynn;
       end
   end
      
   %Update states
   y = ynn;
   t = t + delt;
   %Update output matrix
   if abs(t-outtime(nout+1)) < difftime
      nout = nout+1;
      T2(nout) = t;
      Y2(nout,1) = y(1); Y2(nout,2) = y(2);
   end
end
toc

%% Runge Kutta method
% Approach based on common fourth-order Runge-Kutta method from Wikipedia
disp ('Runge Kutta');

tic
y = yini; %initial states;
t = 0;
%store initial values in output matricese
T3(1) = 0;
Y3(1,1:2) = yini';
nout = 1;

while abs(t-tend) > difftime
   %Calc Rate in order to estimate max timestep
   dy = LotkaVolterraTHe(t,y,Par);
   
   %Check delt
   dttest = abs(0.1*y./dy.*(dy<0))+(dy>=0);
   dtout = outtime(nout+1)-t;
   delt = min([dttest(:)' deltmax dtout]);
   
   %Calc Rates (k1 to k4 for RK4)
   %dy1 = LotkaVolterraTHe(t, y, Par);
   k1 = delt.*dy;
   k2 = delt.*LotkaVolterraTHe(t+delt/2, y+k1/2, Par);
   k3 = delt.*LotkaVolterraTHe(t+delt/2, y+k2/2, Par);
   k4 = delt.*LotkaVolterraTHe(t+delt, y+k3, Par);
      
   %Update states
   y = y + (k1+2*k2+2*k3+k4)/6;
   t = t + delt;
   %Update output matrix
   if abs(t-outtime(nout+1)) < difftime
      nout = nout+1;
      T3(nout) = t;
      Y3(nout,1) = y(1); Y3(nout,2) = y(2);
   end
end
toc

%% Built in ODE Solver
%Using built in ODE solver
options = odeset('RelTol',1e-6,'AbsTol',1e-6);%,'OutputFcn','odeplot');

disp ('Built in ODE solver (ode45)');

tic
[T4,Y4] = ode45(@(t,y) LotkaVolterraTHe(t,y,Par),outtime,yini,options);
toc

%% Built in ODE Solver Implicit!
%Using built in ODE solver ode15i
disp('Built in ODE solver implicit (ode15i)')
options = odeset('RelTol',1e-4,'AbsTol',1e-4);%,'OutputFcn','odeplot');
tic
t0 = 0;
y0 = yini;
yp0 = LotkaVolterraTHe(t0,y0,Par);
%[y0,yp0] = decic(@(t,y,dy) LotkaVolterraImpTHe(t,y,dy,Par),...
%    t0,y0,0,yp0,0);

[T5,Y5] = ode15i(@(t,y,dy) LotkaVolterraImpTHe(t,y,dy,Par),outtime,y0,...
    yp0,options);
toc


figure(1)
clf
plot(T1,Y1,'b.-')
hold on;
plot(T2,Y2,'rx--');
plot(T3,Y3,'k+:');
plot(T4,Y4,'mo-');
plot(T5,Y5,'cp-.');
xlabel('Time')
ylabel('Nr of Animals')
legend({'E_R','E_F', 'P&C_R','P&C_F','RK_R','RK_F','ODE45_R','ODE45_F','ODE15i_R','ODE15i_F'})

figure(2)
clf
plot(Y1(:,1),Y1(:,2),'b.')
hold on
plot(Y2(:,1),Y2(:,2),'ro')
plot(Y3(:,1),Y3(:,2),'g+')
plot(Y4(:,1),Y4(:,2),'m^-')
plot(Y5(:,1),Y5(:,2),'c*:')
xlabel('Rabbits');
ylabel('Foxes');
legend({'Explicit', 'P&C', 'RK','ODE45','ODE15i'})


