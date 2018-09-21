%NL Standing wave (2nd order)

h=0.1; %Step size
N=100; %# of total time steps
x=zeros(N,4);%initialze array
x(1,:)=[1  0]; %initial conditions
ts=0:h:(N-1)*h;

w = @(t) (1+cos(t));
w = @(t) 1;


%x1 will be x and x2 will be dx/dt

f = @(a,b,t) ([b -w(t)*a]); %Slope function for this case

for n=2:N
    
    x1 = x(n-1,1);
    x2 = x(n-1,2);
    th = h*(n-1);
    
    k1 = h.*f(x1,x2,th);
    k2 = h.*f(x1+k1(1,1)/2,x2+k1(1,2)/2,th+h/2);
    k3 = h.*f(x1+k2(1,1)/2,x2+k2(1,2)/2,th+h/2);
    k4 = h.*f(x1+k3(1,1),x2+k3(1,2),th+h);
    
    x(n,:) = x(n-1,:) + 1/6.*(k1+2.*k2+2.*k3+k4);  
    
end