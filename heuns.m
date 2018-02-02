function [t,Vout] = heuns(func,t0,i0,tf,h,Vin,R) 
% solve the ODE i'=func(t,i,Vin) on the interval t0,tf, in steps of h with
% initial values t0,i0
% then convert to Vout by using relationship between i and Vout
N=round((tf-t0)/h);                     % determine size of arrays
ia=zeros(1,N);ta=zeros(1,N);            % set up arrays, initialize values to zero
Vouta=zeros(1,N);                       % set up array for Vout, initialize values to zero
ta(1)=t0;ia(1)=i0;                      % initialize arrays

for j=1:N-1                             % loop for N-1 steps to find Vout(1)...Vout(N-1)
    tt=ta(j);it=ia(j);                  % temporary name for t and i
    Vouta(j)=Vin(j)-R*ia(j);            % find Vout from Vin and i by given equation
    grad1=feval(func,tt,it,Vin(j));     % gradient at t
    ip=it+h*grad1;                      % calculate i-predictor
    grad2=feval(func,tt+h,ip,Vin(j));   % calculate gradient at t+h
    ia(j+1)=it+h*(grad1+grad2)/2;       % next value of i calculated from previous values of t,i
    ta(j+1)=tt+h;                       % increase t by step-size
end

t=ta;Vout=Vouta;                        % return arrays t,Vout