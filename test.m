t0=0;                                  % set initial value of t=0
i0=0;                                  % set initial condition of i=0
tf=0.00013;                               % set final value of t
R=0.5;                                 % set constant value R=0.5
L=0.0015;                              % set constant value L=0.0015
Vin0=5;                                % set initial value of Vin=5
T=0.00013;
func=@(t,i,Vin) Vin/L-R*i/L;
k=R/L;
f=2*pi/T;
g=5/L;

h0=0.0000001;hs=0.0000001;hf=0.00001;
h=h0:hs:hf;
Nh=((hf-h0)/hs)+1;

errora=zeros(1,Nh);

for i=1:Nh
    
    N=round((tf-t0)/h(i));
    Vin=zeros(1,N);ts=zeros(1,N);
    ts(1)=t0;
    
    for j=1:N
      ts1=ts(j);
      Vin(j)=Vin0*cos(2*pi*ts1/T);
      ts(j+1)=ts1+h(i);
    end
    [ta,Vout]=heuns(func,t0,i0,tf,h(i),Vin,R);
    
    ie=((g*f)/(f^2+k^2))*(sin(f*ta)+(k/f)*cos(f*ta))-(g*k)/(f^2+k^2);
    exact=Vin-R*ie;
    error=abs(exact-Vout);
    errora(i)=max(error);
end

loglog(h,errora);
    