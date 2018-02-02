function [t,Vout] = midpoint(func,t0,i0,tf,h,Vin,R)
% solve the ODE i'=func(t,i,Vin) on the interval t0,tf, in steps of h with
% initial values t0,i0
% then convert to Vout by using relationship between i and Vout
N=round((tf-t0)/h);                     % determine size of arrays
ia=zeros(1,N);ta=zeros(1,N);            % set up arrays, initialize values to zero
Vouta=zeros(1,N);                       % set up array for Vout, initialize values to zero
ta(1)=t0;ia(1)=i0;                      % initialize arrays

for j=1:N-1                             % loop for N-1 steps to find Vout(1)...Vout(N-1)
    k1=feval(func,ta(j),ia(j),Vin(j));  % gradient at t
    ia(j+1)=ia(j)+h*feval(func,ta(j)+1/2*h,ia(j)+1/2*k1*h,Vin(j));
    % next value of i calculated from previous values of t,i
    ta(j+1)=ta(j)+h;                    % increase t by step-size
    Vouta(j)=Vin(j)-R*ia(j);            % find Vout from Vin and i by given equation
end

t=ta;Vout=Vouta;                        % return arrays t,Vout