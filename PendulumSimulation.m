% Lars Osterberg
%
% For all my code please run it section by section
%
%% Euler's method results compared with runge Kutta

% Parameters
h=0.001; % Step size
t_final=20; % Final time
n_steps=t_final/h; % Number of time steps
m=0.23; % mass in kg
g=9.81; % gravitational constant in m/s^2
L=0.1549; % length in meters from pivot to center of mass
r=0.0127/2; % radius in meters
mu=0.11; % coefficient of friction
I=0.00896; % moment of inertia calculated using rigid body
% Initial conditions
y1=-pi/4; % Initial angle 
y2=0; % Initial velocity
% Arrays to store results
teuler=0:h:t_final; % Time vector for euler's method
y1_vals=zeros(1, n_steps+1); % y1 values
y2_vals=zeros(1, n_steps+1); % y2 values
% Initial values
y1_vals(1)=y1;
y2_vals(1)=y2;

% Euler's Method
for i=1:n_steps
    % Current values
    y1_curr=y1_vals(i);
    y2_curr=y2_vals(i);
    % Compute derivatives
    dy1_dt=y2_curr;
    dy2_dt=-(L*g*m*sin(y1_curr)/I)-(mu*r*(g*m*cos(y1_curr)+L*m*y2_curr^2)/I)*sign(y2_curr); %equation of motion   
    % Euler update
    y1_vals(i+1)=y1_curr+h*dy1_dt;
    y2_vals(i+1)=y2_curr+h*dy2_dt;
end

% now for Runge Kutta
[tRK,yRK] = ode45(@thetaRK,[0 20],[-pi/4; 0]); % performing the numerical simulation for Runge Kutta
% Plotting results
figure;
hold on
plot(teuler, y1_vals*180/pi); %plotting the results
plot(tRK,yRK(:,1)*180/pi) %plotting the results
hold off
title('Comparing Euler with Runge-Kutta for 1/2" bushing');
xlabel('Time (seconds)');
ylabel('Angle (degrees)');
legend('Euler','Runge Kutta')

function dydt = thetaRK(t,y)
    m=0.23;
    g=9.81;
    L=0.1549;
    r=0.0127/2;
    mu=0.11;
    I=0.00896;
    dydt=[y(2); -(L*g*m*sin(y(1))/I)-(mu*r/I)*(L*m*y(2)^2+g*m*cos(y(1)))*sign(y(2))]; %equation of motion 
end

%% Part A results verification
[tA,yA] = ode45(@thetaA,[0 20],[-pi/16; 0]); %performing the numerical simulation using a small angle (pi/16 is small enough)
plot(tA,yA(:,1)) %plotting the results
%by observing this plot, we see the period is about 0.804 seconds

function dydt = thetaA(t,y)
    m=0.23;
    g=9.81;
    L=0.1549;
    r=0.0127;
    mu=0; % no friction
    I=m*L^2; % point mass assumption
    dydt=[y(2); -(L*g*m*sin(y(1))/I)-(mu*r/I)*(L*m*y(2)^2+g*m*cos(y(1)))*sign(y(2))]; %equation of motion  
end
%% Part B results verification
[tB,yB] = ode45(@thetaB,[0 20],[-pi/16; 0]); %performing the numerical simulation using a small angle (pi/16 is small enough)
plot(tB,yB(:,1)) %plotting the results
%by observing this plot, we see the period is about 1.007 seconds

function dydt = thetaB(t,y)
    m=0.23;
    g=9.81;
    L=0.1549;
    r=0.0127;
    mu=0; % no friction
    I=0.00896; %moment of inertia calculated by hand using rigid body
    dydt=[y(2); -(L*g*m*sin(y(1))/I)-(mu*r/I)*(L*m*y(2)^2+g*m*cos(y(1)))*sign(y(2))]; %equation of motion   
end
%% Part C results verification
[tC,yC] = ode45(@thetaC,[0 20],[-pi/4; 0]); %performing the numerical simulation for 45 degree starting angle 
plot(tC,yC(:,1)) %plotting the results
%by observing this plot, we see the period is about 1.035 seconds

function dydt = thetaC(t,y)
    m=0.23;
    g=9.81;
    L=0.1549;
    r=0.0127;
    mu=0; %no friction
    I=0.00896; % moment of inertia calculated using rigid body
    dydt=[y(2); -(L*g*m*sin(y(1))/I)-(mu*r/I)*(L*m*y(2)^2+g*m*cos(y(1)))*sign(y(2))]; %equation of motion   
end
%% Part D results verification
options = odeset('Events', @myEventsFcn); %Here I am using an event function to find total oscillation time
[thalf,yhalf,totalhalf] = ode45(@thetaDhalf,[0 50],[-pi/4; 0],options); %performing the numerical simulation 1/2" bushing
[tquarter,yquarter,totalquarter] = ode45(@thetaDquarter,[0 90],[-pi/4; 0],options); %performing the numerical simulation for 1/4"
[teighth,yeighth,totaleighth] = ode45(@thetaDeighth,[0 160],[-pi/4; 0],options); %performing the numerical simulation 1/8" bushing
%the total time for each bushing are the variables totalhalf,
%totalquarter, and totaleighth

function dydt = thetaDhalf(t,y) %1/2" bushing
    m=0.23;
    g=9.81;
    L=0.1549;
    r=0.0127/2;
    mu=0.11;
    I=0.00896;
    dydt=[y(2); -(L*g*m*sin(y(1))/I)-(mu*r/I)*(L*m*y(2)^2+g*m*cos(y(1)))*sign(y(2))]; %equation of motion 
end
function dydt = thetaDquarter(t,y) %1/4" bushing
    m=0.23;
    g=9.81;
    L=0.1549;
    r=0.0127/4;
    mu=0.11;
    I=0.00896;
    dydt=[y(2); -(L*g*m*sin(y(1))/I)-(mu*r/I)*(L*m*y(2)^2+g*m*cos(y(1)))*sign(y(2))]; %equation of motion 
end
function dydt = thetaDeighth(t,y) %1/8" bushing
    m=0.23;
    g=9.81;
    L=0.1549;
    r=0.0127/8;
    mu=0.11;
    I=0.00896;
    dydt=[y(2); -(L*g*m*sin(y(1))/I)-(mu*r/I)*(L*m*y(2)^2+g*m*cos(y(1)))*sign(y(2))]; %equation of motion 
end
function [position,isterminal,direction] = myEventsFcn(t,y)
  position = abs(y(1))+abs(y(2))-0.01; % pendulum is considered stopped when this equals 0, or abs(y(1))+abs(y(2))<0.01
  isterminal = 1;  % Halt integration at this point
  direction = 0;   % The zero can be approached from either direction
end
%% Plotting for 1/4" bushing
options = odeset('Events', @myEventsFcn2); %Here I am using an event function to find total oscillation time
[t,y,te,ye,ie] = ode45(@thetaquarterbushing,[0 90],[-pi/4; 0],options); %performing the numerical simulation

figure
plot(t,180*y(:,1)/pi) %plotting the results for the angle
xline(te(1,1)) %plotting where it stops using my events function, te is the end time
xlabel('Time (seconds)')
ylabel('Angle (degrees)')
title(['Simulation of oscillation angle for 1/4" bushing pendulum, stop time at ', num2str(te(1,1)), ' seconds'])

figure
plot(t,180*y(:,2)/pi) %plotting the results for the angular velocity
xline(te(1,1)) %plotting where it stops using my events function, te is the end time
xlabel('Time (seconds)')
ylabel('Angular Velocity (degrees/second)')
title(['Simulation of angular velocity for 1/4" bushing pendulum, stop time at ', num2str(te(1,1)), ' seconds'])

figure
m=0.23;
g=9.81;
L=0.1549;
PE=m*g*L*(1-cos(y(:,1))); %potential energy
plot(t,PE) %plotting the results for the potential energy
xline(te(1,1)) %plotting where it stops using my events function, te is the end time
xlabel('Time (seconds)')
ylabel('Potential Energy (joules)')
title(['Simulation of potential energy for 1/4" bushing pendulum, stop time at ', num2str(te(1,1)), ' seconds'])

figure
I=0.00896;
KE=0.5*I*(y(:,2)).^2; %kinetic energy
plot(t,KE) %plotting the results for the kinetic energy
xline(te(1,1)) %plotting where it stops using my events function, te is the end time
xlabel('Time (seconds)')
ylabel('Kinetic Energy (joules)')
title(['Simulation of kinetic energy for 1/4" bushing pendulum, stop time at ', num2str(te(1,1)), ' seconds'])

figure
TE=PE+KE; %total energy
plot(t,TE) %plotting the results for the total energy
xline(te(1,1)) %plotting where it stops using my events function, te is the end time
xlabel('Time (seconds)')
ylabel('Total Energy (joules)')
title(['Simulation of total energy for 1/4" bushing pendulum, stop time at ', num2str(te(1,1)), ' seconds'])

function dydt = thetaquarterbushing(t,y) %bushing with r=1/4"
    m=0.23;
    g=9.81;
    L=0.1549;
    r=0.0127/4;
    mu=0.11;
    I=0.00896;
    dydt=[y(2); -(L*g*m*sin(y(1))/I)-(mu*r/I)*(L*m*y(2)^2+g*m*cos(y(1)))*sign(y(2))]; %equation of motion  
end
function [position,isterminal,direction] = myEventsFcn2(t,y)
  position = abs(y(1))+abs(y(2))-0.01; % pendulum is considered stopped when this equals 0, or abs(y(1))+abs(y(2))<0.01
  isterminal = 0;  % Do not halt integration at this point
  direction = 0;   % The zero can be approached from either direction
end