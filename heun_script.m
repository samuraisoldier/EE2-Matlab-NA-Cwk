% Case 1
t0=0;                                  % set initial value of t=0
i0=0;                                  % set initial condition of i=0
tf=0.03;                               % set final value of t
h=0.0007;                              % set step-size
N=round((tf-t0)/h);                    % set size of arrays

R=0.5;                                 % set constant value R=0.5
L=0.0015;                              % set constant value L=0.0015
Vin0=5;                                % set initial value of Vin=5

func=@(t,i,Vin) Vin/L-R*i/L;           % Li'+Ri=Vin -> i'=Vin/L-Ri/L

Vin=zeros(1,N);                        % set up arrays for Vin
for j=1:N-1                            % loop for N steps
    Vin(j)=Vin0;                       % set all elements in Vin=5
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);% call heuns.m

figure                                 % new figure
subplot(2,1,1);                        % plotting the first graph on the figure
plot(ta,Vout);                         % plot approximation Vout/t
grid on                                % set up major grid lines in graph
grid minor                             % set up minor grid lines in graph
title('Vout/t for RL circuit with Vin=5V') % set up title for the graphs
xlabel('0 < t < 0.03')                 % label the x-axis
ylabel('Vout')                         % label the y-axis
hold on;                               % wait for next graph to be plot on the same axes
exact=Vin0-R*((Vin0/R)*(1-exp(-(R/L)*ta))); % exact solution of ODE
plot(ta,exact);                        % plot exact on the same axes with approximation
legend('approximation','exact')        % set up legend for approximation and exact
subplot(2,1,2);                        % plot second graph on the figure
plot(ta,exact-Vout);                   % plot error as a function of t
grid on
grid minor
title('Error in Vout for RL circuit with Vin=5V')
xlabel('0 < t < 0.03')
ylabel('Error')
max(abs(exact-Vout));                  % calulate max value of error over range of t

% Case 2a1
tf=0.01;                               % set up final t=0.01
h=0.0001;                              % set up step-size
N=round((tf-t0)/h);                    % new array size

Vin0=4;                                % new initial value of Vin=4
tau1=0.00014^2;                        % constant value tau=0.00014^2

Vin=zeros(1,N);ts=zeros(1,N);          % set up arrays
ts(1)=t0;                              % first element in temporary array for 
                                       % t equals initial value of t
for j=1:N-1                            % loop for N steps
    ts1=ts(j);                         % temporary name for t
    Vin(j)=Vin0*exp(-ts(j)^2/tau1);    % set up Vin value at given t
    ts(j+1)=ts1+h;                     % increase t by step-size
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);% call heuns.m

figure;
subplot(2,2,1);
plot(ta,Vout,'r*')
grid on
grid minor
title('Vout/t for RL circuit with \itVin=4e^{-t^2/\tau}')
xlabel('0 < t < 0.01')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b')
legend('Vout','Vin')

% Case 2a2
tf=0.05;                               % new final value for t
h=0.0005;                              % set up step-size
N=round((tf-t0)/h);

Vin0=4;
tau2=0.00014;                          % new tau value

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N-1
    ts1=ts(j);
    Vin(j)=Vin0*exp(-ts(j)^2/tau2);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

subplot(2,2,3);
plot(ta,Vout,'r*')
grid on
grid minor
title('Vout/t for RL circuit with \itVin=4e^{-t^2/\tau}')
xlabel('0 < t < 0.05')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b')
legend('Vout','Vin')

% Case 2b1
tf=0.007;
h=0.0001;
N=round((tf-t0)/h);

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N-1
    ts1=ts(j);
    Vin(j)=Vin0*exp(-ts(j)/tau1);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

subplot(2,2,2);
plot(ta,Vout,'r*')
grid on
grid minor
title('Vout/t for RL circuit with \itVin=4e^{-t/\tau}')
xlabel('0 < t < 0.007')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b')
legend('Vout','Vin')

% Case 2b2
tf=0.01;
h=0.0001;
N=round((tf-t0)/h);

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N-1
    ts1=ts(j);
    Vin(j)=Vin0*exp(-ts(j)/tau2);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

subplot(2,2,4);
plot(ta,Vout,'r*');
grid on
grid minor
title('Vout/t for RL circuit with \itVin=4e^{-t/\tau}')
xlabel('0 < t < 0.01')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b');
legend('Vout','Vin')

% Case 3a1
tf=0.00014;
h=0.000001;
N=round((tf-t0)/h);

Vin0=5;
T1=0.00014;

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N
    ts1=ts(j);
    Vin(j)=Vin0*sin(2*pi*ts1/T1);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

figure;
subplot(2,2,1);
plot(ta,Vout,'r*');
grid on
grid minor
title('Vout/t for RL circuit with \itVin as sinewave with T=140\mus')
xlabel('0 < t < 140\mus')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b');
legend('Vout','Vin')

% Case 3a2
tf=0.000025;
h=0.0000005;
N=round((tf-t0)/h);

T2=0.000025;

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N
    ts1=ts(j);
    Vin(j)=Vin0*sin(2*pi*ts1/T2);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

subplot(2,2,2);
plot(ta,Vout,'r*');
grid on
grid minor
title('Vout/t for RL circuit with \itVin as sinewave with T=25\mus')
xlabel('0 < t < 25\mus')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b');
legend('Vout','Vin')

% Case 3a3
tf=0.00055;
h=0.000005;
N=round((tf-t0)/h);

T3=0.00055;

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N
    ts1=ts(j);
    Vin(j)=Vin0*sin(2*pi*ts1/T3);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

subplot(2,2,3);
plot(ta,Vout,'r*');
grid on
grid minor
title('Vout/t for RL circuit with \itVin as sinewave with T=550\mus')
xlabel('0 < t < 550\mus')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b');
legend('Vout','Vin')

% Case 3a4
tf=0.0012;
h=0.00001;
N=round((tf-t0)/h);

T4=0.0012;

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N
    ts1=ts(j);
    Vin(j)=Vin0*sin(2*pi*ts1/T4);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

subplot(2,2,4);
plot(ta,Vout,'r*');
grid on
grid minor
title('Vout/t for RL circuit with \itVin as sinewave with T=1200\mus')
xlabel('0 < t < 1200\mus')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b');
legend('Vout','Vin')

% Case 3b1
tf=0.00014;
h=0.000001;
N=round((tf-t0)/h);

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N
    ts1=ts(j);
    Vin(j)=Vin0*square(2*pi*ts1/T1);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

figure
subplot(2,2,1);
plot(ta,Vout,'r*');
grid on
grid minor
title('Vout/t for RL circuit with \itVin as squarewave with T=140\mus')
xlabel('0 < t < 140\mus')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b');
legend('Vout','Vin')

% Case 3b2
tf=0.000025;
h=0.0000005;
N=round((tf-t0)/h);

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N
    ts1=ts(j);
    Vin(j)=Vin0*square(2*pi*ts1/T2);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

subplot(2,2,2);
plot(ta,Vout,'r*');
grid on
grid minor
title('Vout/t for RL circuit with \itVin as squarewave with T=25\mus')
xlabel('0 < t < 25\mus')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b');
legend('Vout','Vin')

% Case 3b3
tf=0.00055;
h=0.000005;
N=round((tf-t0)/h);

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N
    ts1=ts(j);
    Vin(j)=Vin0*square(2*pi*ts1/T3);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

subplot(2,2,3);
plot(ta,Vout,'r*');
grid on
grid minor
title('Vout/t for RL circuit with \itVin as squarewave with T=550\mus')
xlabel('0 < t < 550\mus')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b');
legend('Vout','Vin')

% Case 3b4
tf=0.0012;
h=0.00001;
N=round((tf-t0)/h);

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N
    ts1=ts(j);
    Vin(j)=Vin0*square(2*pi*ts1/T4);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

subplot(2,2,4);
plot(ta,Vout,'r*');
grid on
grid minor
title('Vout/t for RL circuit with \itVin as squarewave with T=1200\mus')
xlabel('0 < t < 1200\mus')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b');
legend('Vout','Vin')

%Case 3c1
tf=0.00028;
h=0.000001;
N=round((tf-t0)/h);

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N
    ts1=ts(j);
    Vin(j)=Vin0*sawtooth(2*pi*ts1/T1);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

figure
subplot(2,2,1);
plot(ta,Vout,'r*');
grid on
grid minor
title('Vout/t for RL circuit with \itVin as sawtooth with T=140\mus')
xlabel('0 < t < 280\mus')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b');
legend('Vout','Vin')

% Case 3c2
tf=0.00005;
h=0.0000005;
N=round((tf-t0)/h);

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N
    ts1=ts(j);
    Vin(j)=Vin0*sawtooth(2*pi*ts1/T2);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

subplot(2,2,2);
plot(ta,Vout,'r*');
grid on
grid minor
title('Vout/t for RL circuit with \itVin as sawtooth with T=25\mus')
xlabel('0 < t < 50\mus')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b');
legend('Vout','Vin')

% Case 3c3
tf=0.00055*2;
h=0.000005;
N=round((tf-t0)/h);

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N
    ts1=ts(j);
    Vin(j)=Vin0*sawtooth(2*pi*ts1/T3);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

subplot(2,2,3);
plot(ta,Vout,'r*');
grid on
grid minor
title('Vout/t for RL circuit with \itVin as sawtooth with T=550\mus')
xlabel('0 < t < 1100\mus')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b');
legend('Vout','Vin')

% Case 3c4
tf=0.0012*2;
h=0.00001;
N=round((tf-t0)/h);

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N
    ts1=ts(j);
    Vin(j)=Vin0*sawtooth(2*pi*ts1/T4);
    ts(j+1)=ts1+h;
end

[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);

subplot(2,2,4);
plot(ta,Vout,'r*');
grid on
grid minor
title('Vout/t for RL circuit with \itVin as sawtooth with T=1200\mus')
xlabel('0 < t < 2400\mus')
ylabel('Voltage')
hold on;
plot(ta,Vin,'b');
legend('Vout','Vin')