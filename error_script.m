% Heun's Method
t0=0;                                  % set initial value of t=0
i0=0;                                  % set initial condition of i=0
tf=0.00013;                               % set final value of t
h=0.000001;                              % set step-size
N=round((tf-t0)/h);                    % set size of arrays

R=0.5;                                 % set constant value R=0.5
L=0.0015;                              % set constant value L=0.0015
Vin0=5;                                % set initial value of Vin=5
T=0.00013;

func=@(t,i,Vin) Vin/L-R*i/L;           % Li'+Ri=Vin -> i'=Vin/L-Ri/L

Vin=zeros(1,N);ts=zeros(1,N);
ts(1)=t0;
for j=1:N
    ts1=ts(j);
    Vin(j)=Vin0*cos(2*pi*ts1/T);
    ts(j+1)=ts1+h;
end


[ta,Vout]=heuns(func,t0,i0,tf,h,Vin,R);
k=R/L;
f=2*pi/T;
g=5/L;
ie=((g*f)/(f^2+k^2))*(sin(f*ta)+(k/f)*cos(f*ta))-(g*k)/(f^2+k^2);
exact=Vin-R*ie;
figure
subplot(3,3,1)
plot(ta,Vout,'r*');
grid on                                % set up major grid lines in graph
grid minor                             % set up minor grid lines in graph
title('Vout/t for RL circuit with Vin=consine with period 130\mus') % set up title for the graphs
xlabel('0 < t < 130\mus')              % label the x-axis
ylabel('Vout')                         % label the y-axis
hold on
plot(ta,exact,'b');
legend('approximation','exact')        % set up legend for approximation and exact
error=abs(exact-Vout);
subplot(3,3,4)
plot(ta,error);                   % plot error as a function of t
grid on
grid minor
title('Error in Vout for RL circuit with Vin=cosin with period 130\mus')
xlabel('0 < t < 130\mus')
ylabel('Error')
indexmax = find(max(error) == error);
tmax = ta(indexmax);
errormax = error(indexmax);
strmax = ['Maximum = ', num2str(errormax)];
text(tmax,errormax,strmax,'HorizontalAlignment','right');
hold on
plot(tmax,errormax,'rX')



