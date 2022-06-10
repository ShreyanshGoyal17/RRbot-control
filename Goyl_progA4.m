clc;clear all
syms theta1 theta2 theta_dot1 theta_dot2 theta_doubledot1 theta_doubledot2 g t 'real'
syms m1 m2 l1 l2 r1 r2 I1 I2 'real'


% the equations of motion can be transformed into the standard Manipulator form:
% M (q) ̈q + C(q,  ̇q)  ̇q + g(q) = τ


gq=[-g*m2*(r2*sin(theta1 + theta2) + l1*sin(theta1)) - g*m1*r1*sin(theta1);-g*m2*r2*sin(theta1 + theta2)];

M1=  [m1*r1^2 + I1 + I2 + (m2*(2*l1^2 + 4*cos(theta2)*l1*r2 + 2*r2^2))/2, I2 + (m2*(2*r2^2 + 2*l1*cos(theta2)*r2))/2];
 
M2= [I2 + (m2*(2*r2^2 + 2*l1*cos(theta2)*r2))/2, m2*r2^2 + I2];
 
Mq= [M1;M2];

Cq_qdot= [-(m2*theta_dot2*(2*l1*r2*sin(theta2)*(theta_dot1 + theta_dot2) + 2*l1*r2*theta_dot1*sin(theta2)))/2;...
    l1*m2*r2*theta_dot1*sin(theta2)*(theta_dot1 + theta_dot2) - l1*m2*r2*theta_dot1*theta_dot2*sin(theta2)];
 
u= Mq*[theta_doubledot1;theta_doubledot2] + Cq_qdot+ gq;


%the given nominal values

m_hat=0.75; m2_hat=0.75; I1_hat=0.063; I2_hat= 0.063;

A= [0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
B=[0 0;0 0;1 0;0 1];

lamda= [-3 -5 -5 -7];

K=place(A,B,lamda);

Asys=[0 0 1 0;0 0 0 1;-K];

%calculate P
Q= diag([10 10 10 10]);
P=lyap(Asys',Q);

%running the ode45 on the function ode_system defined below
T=10;
y0=[deg2rad(200),deg2rad(125),0,0];
[t,y]=ode45(@ode_system,[0,T],y0);

rho= 2;


% getting control inputs given by the ode function below

Tau=[];
for index=1:length(t)
    time=t(index);
    x=y(index,:).';
    [~,tau]=ode_system(time,x);
    Tau=[Tau,tau];
end



%the given desired trajectories
 
q1_d= ((pi*t.^3)/500) + (-3*pi*t.^2)/100 + pi;
q1dot_d= ((3*pi*t.^2)/500) + 6*(-pi*t)/100;
q1ddot_d= ((6*pi*t)/500) - 6*pi/100 ;

q2_d= ((pi*t.^3)/1000) -(3*pi*t.^2)/200 + pi/2;
q2dot_d= ((3*pi*t.^2)/1000) + 6*(-pi*t)/200;
q2ddot_d= ((6*pi*t)/1000) - 6*pi/200;
qd= [q1_d ;q2_d ;q1dot_d; q2dot_d];
x= [theta1-q1_d ;theta2-q2_d ;theta_dot1-q1dot_d; theta_dot2-q2dot_d]; %states to control


figure;

subplot(2,2,1)
hold on
plot(t,y(:,1));
xlabel('Time')
ylabel('theta1')
plot(t,q1_d)

subplot(2,2,2)
hold on
plot(t,y(:,2));
xlabel('Time')
ylabel('theta2')
plot(t,q2_d)

subplot(2,2,3);
hold on
plot(t,y(:,3));
xlabel('Time')
ylabel('theta1 dot')
plot(t,q1dot_d)

subplot(2,2,4);
hold on
plot(t,y(:,4));
xlabel('Time')
ylabel('theta2 dot')
plot(t,q2dot_d)

figure;
 
subplot(2,2,1)
plot(t,Tau(1,:))
xlabel('Time')
ylabel('Tau1')
subplot(2,2,2)
plot(t,Tau(2,:))
xlabel('Time')
ylabel('Tau2')


%this function outputs states x and torque tau at each time instance with the help of ode45

function [dx,tau]=ode_system(t,x)

m1=1; m2=1; l1=1; l2=1; r1=0.45; r2= 0.45; I1=0.084; I2=0.084; g=9.81;
m1_hat=0.75; m2_hat=0.75; I1_hat=0.063; I2_hat= 0.063;

dx= zeros(4,1);
x=num2cell(x);
[theta1,theta2,theta_dot1,theta_dot2]=deal(x{:});

gq=[- g*m2_hat*(r2*sin(theta1 + theta2) + l1*sin(theta1)) - g*m1_hat*r1*sin(theta1);-g*m2_hat*r2*sin(theta1 + theta2)];

M1=  [m1_hat*r1^2 + I1_hat + I2_hat + (m2_hat*(2*l1^2 + 4*cos(theta2)*l1*r2 + 2*r2^2))/2, I2_hat + (m2_hat*(2*r2^2 + 2*l1*cos(theta2)*r2))/2];
 
M2= [I2_hat + (m2_hat*(2*r2^2 + 2*l1*cos(theta2)*r2))/2, m2_hat*r2^2 + I2_hat];
 
Mq= [M1;M2];

Cq_qdot= [-(m2_hat*theta_dot2*(2*l1*r2*sin(theta2)*(theta_dot1 + theta_dot2) + 2*l1*r2*theta_dot1*sin(theta2)))/2;...
    l1*m2_hat*r2*theta_dot1*sin(theta2)*(theta_dot1 + theta_dot2) - l1*m2_hat*r2*theta_dot1*theta_dot2*sin(theta2)];
 
 
%desired polynomial path for joint 1 and 2

q1_d= ((pi*t.^3)/500) + (-3*pi*t.^2)/100 + pi;
q1dot_d= ((3*pi*t.^2)/500) + 6*(-pi*t)/100;
q1ddot_d= ((6*pi*t)/500) - 6*pi/100 ;

q2_d= ((pi*t.^3)/1000) -(3*pi*t.^2)/200 + pi/2;
q2dot_d= ((3*pi*t.^2)/1000) + 6*(-pi*t)/200;
q2ddot_d= ((6*pi*t)/1000) - 6*pi/200;
qd= [q1_d ;q2_d ;q1dot_d; q2dot_d];

A= [0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
B=[0 0;0 0;1 0;0 1];

K=[35 0 12 0;0 15 0 8];
Asys=[0 0 1 0;0 0 0 1;-K];

%calculate P
Q= diag([10 10 10 10]);
P=lyap(Asys',Q);

rho= 2;
x= [theta1-q1_d ;theta2-q2_d ;theta_dot1-q1dot_d; theta_dot2-q2dot_d];
phi=0.09;

%comment the below section to run without the boundary condition
%boundary condition to avoid chattering
if phi>0
    if norm(B'*P*x)>phi
        vr=-(rho*B'*P*x)/norm(B'*P*x);
    else
        vr=-(rho*B'*P*x)/phi;
    end

else
  if norm(B'*P*x)~=0
    vr= -(rho*B'*P*x)/norm(B'*P*x);
  else
    vr=0;
  end
end


v= -K*x + [q1ddot_d;q2ddot_d] + vr;

tau= Mq*v+ Cq_qdot+ gq;

u1=tau(1);
u2=tau(2);
dx(1)=theta_dot1;
dx(2)= theta_dot2;
dx(3)= (I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 - l1*m2^2*r2^3*theta_dot1^2*sin(theta2) - l1*m2^2*r2^3*theta_dot2^2*sin(theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u2*cos(theta2) - 2*l1*m2^2*r2^3*theta_dot1*theta_dot2*sin(theta2) + l1^2*m2^2*r2^2*theta_dot1^2*cos(theta2)*sin(theta2) + g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) - I2*l1*m2*r2*theta_dot1^2*sin(theta2) - I2*l1*m2*r2*theta_dot2^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) - 2*I2*l1*m2*r2*theta_dot1*theta_dot2*sin(theta2))/(I1*I2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + I1*m2*r2^2 - l1^2*m2^2*r2^2*cos(theta2)^2);
dx(4)=(I1*u2 - I2*u1 + I2*u2 + l1^2*m2*u2 - m2*r2^2*u1 + m2*r2^2*u2 + l1*m2^2*r2^3*theta_dot1^2*sin(theta2) + l1^3*m2^2*r2*theta_dot1^2*sin(theta2) + l1*m2^2*r2^3*theta_dot2^2*sin(theta2) + g*l1^2*m2^2*r2*sin(theta1 + theta2) + I1*g*m2*r2*sin(theta1 + theta2) - g*l1*m2^2*r2^2*sin(theta1) - I2*g*l1*m2*sin(theta1) - I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta_dot1*theta_dot2*sin(theta2) - 2*l1^2*m2^2*r2^2*theta_dot1^2*cos(theta2)*sin(theta2) - l1^2*m2^2*r2^2*theta_dot2^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*m2^2*r2*cos(theta2)*sin(theta1) + I1*l1*m2*r2*theta_dot1^2*sin(theta2) + I2*l1*m2*r2*theta_dot1^2*sin(theta2) + I2*l1*m2*r2*theta_dot2^2*sin(theta2) - g*m1*m2*r1*r2^2*sin(theta1) - 2*l1^2*m2^2*r2^2*theta_dot1*theta_dot2*cos(theta2)*sin(theta2) + 2*I2*l1*m2*r2*theta_dot1*theta_dot2*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(I1*I2 + l1^2*m2^2*r2^2 + I2*l1^2*m2+ I1_hat*m2*r2^2 - l1^2*m2^2*r2^2*cos(theta2)^2);

end




