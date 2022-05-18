
clear; close; clc;
%ROS Setup
rosinit;

j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');
tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);
tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(200), deg2rad(125)];
resp = call(client,req,'Timeout',3);
% tic;
m1=1; m2=1; l1=1; l2=1; r1=0.45; r2= 0.45; I1=0.084; I2=0.084; g=9.81;
m1_hat=0.75; m2_hat=0.75; I1_hat=0.063; I2_hat= 0.063;

A= [0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
B=[0 0;0 0;1 0;0 1];

%K=[12 0 7 0;0 12 0 7];
K=[35 0 12 0;0 15 0 8];
Asys=[0 0 1 0;0 0 0 1;-K];

%calculate P
Q= diag([20 20 20 20]);
P=lyap(Asys',Q);
x= zeros(4,1);
theta1=x(1);
theta2=x(2);
theta_dot1=x(3);
theta_dot2=x(4);

rho= 5;
t = 0;
tau1_w=[];
tau2_w=[];
position=[];
velocity=[];
time=[];

i=1;
tic;
while(t < 10)
t = toc;
%read the joint states
jointData = receive(JointStates);
% inspect the "jointData" variable in MATLAB to get familiar with its structure
% design your state feedback controller in the following
theta1=jointData.Position(1);
theta2=jointData.Position(2);
theta_dot1=jointData.Velocity(1);
theta_dot2=jointData.Velocity(2);

gq=[-g*m2_hat*(r2*sin(theta1 + theta2) + l1*sin(theta1)) - g*m1_hat*r1*sin(theta1);-g*m2_hat*r2*sin(theta1 + theta2)];
%Mq=[I1_hat+I2_hat+m2_hat*l1^2+m2_hat*r2^2-2*l1*m2_hat*cos(theta2) m2_hat*r2^2+ I2_hat-m2_hat*r2*l1*cos(theta2); m2_hat*r2^2+I2_hat- m2_hat*r2*l1*cos(theta2) I2_hat+m2_hat*r2^2];
M1=  [m1_hat*r1^2 + I1_hat + I2_hat + (m2_hat*(2*l1^2 + 4*cos(theta2)*l1*r2 + 2*r2^2))/2, I2_hat + (m2_hat*(2*r2^2 + 2*l1*cos(theta2)*r2))/2];
 
M2= [I2_hat + (m2_hat*(2*r2^2 + 2*l1*cos(theta2)*r2))/2, m2_hat*r2^2 + I2_hat];
 
Mq= [M1;M2];

%Cq_qdot= [l1*m2_hat*r2*sin(theta2)*theta_dot2^2+ 2*l1*m2_hat*r2*theta_dot1*sin(theta2)*theta_dot2;-l1*m2_hat*r2*theta_dot1^2*sin(theta2)];
Cq_qdot= [-(m2_hat*theta_dot2*(2*l1*r2*sin(theta2)*(theta_dot1 + theta_dot2) + 2*l1*r2*theta_dot1*sin(theta2)))/2;...
    l1*m2_hat*r2*theta_dot1*sin(theta2)*(theta_dot1 + theta_dot2) - l1*m2_hat*r2*theta_dot1*theta_dot2*sin(theta2)];

q1_d= ((pi*t.^3)/500) + (-3*pi*t.^2)/100 + pi;
q1dot_d= ((3*pi*t.^2)/500) + 6*(-pi*t)/100;
q1ddot_d= ((6*pi*t)/500) - 6*pi/100 ;

q2_d= ((pi*t.^3)/1000) -(3*pi*t.^2)/200 + pi/2;
q2dot_d= ((3*pi*t.^2)/1000) + 6*(-pi*t)/200;
q2ddot_d= ((6*pi*t)/1000) - 6*pi/200;
qd= [q1_d ;q2_d ;q1dot_d; q2dot_d];


x= [theta1-q1_d ;theta2-q2_d ;theta_dot1-q1dot_d; theta_dot2-q2dot_d];
phi=0.09;
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

tau1.Data= tau(1,:);
tau2.Data= tau(2,:);


% tau1.Data = -K(1,:)*[jointData.Position(1);jointData.Position(2);jointData.Velocity(1);jointData.Velocity(2)];
% tau2.Data = -K(2,:)*[jointData.Position(1);jointData.Position(2);jointData.Velocity(1);jointData.Velocity(2)];

send(j1_effort,tau1);
send(j2_effort,tau2);
% you can sample data here to be plotted at the end

time=[time t];
tau1_w=[tau1_w tau1.Data];
tau2_w=[tau2_w tau2.Data];
position=[position jointData.Position];
velocity=[velocity jointData.Velocity];

end

figure
subplot(2,2,1)
plot(time,tau1_w(1,:))
xlabel('Time')
ylabel('Tau1')
subplot(2,2,2)
plot(time,tau2_w(1,:))
xlabel('Time')
ylabel('Tau2')

figure
subplot(2,2,1)
plot(time,position(1,:))
xlabel('Time')
ylabel('Theta1')
subplot(2,2,2)
plot(time,position(2,:))
xlabel('Time')
ylabel('Theta2')
subplot(2,2,3)
plot(time,velocity(1,:))
xlabel('Time')
ylabel('Theta1-dot')
subplot(2,2,4)
plot(time,velocity(2,:))
xlabel('Time')
ylabel('Theta2-dot')


tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);

%disconnect from roscore
rosshutdown;


