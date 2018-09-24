%NL Standing wave (3rd order)
clear all
g = 1;
L = 1;
j = 1;


dt=0.01; %Step size
N=1000; %# of total time steps
x=zeros(N,4);%initialze array
x(1,:)=[0.01 0.05 0.05 0.01];
Aj = x(1,1)/2 - 1i*x(1,2)*k(j)*T(j)/(2*omega(j));
ts=0:dt:(N-1)*dt;
H=zeros(1,N);
H(1,1) = hamiltonian(x(1,:),j);

%RK4 method

for n=2:N
    
    x_old = x(n-1,:);
    
    k1 = dt.*f(x_old,j);
    k2 = dt.*f(x_old+k1./2,j);
    k3 = dt.*f(x_old+k2./2,j);
    k4 = dt.*f(x_old+k3,j);
    
    x(n,:) = x_old + 1/6.*(k1+2.*k2+2.*k3+k4); 
    H(1,n) = hamiltonian(x(n,:),j);
    
end

%Plot harmonics and phase planes
figure()
subplot(2,2,1), plot(ts,x(:,1))
hold on
plot(ts,x(:,2))
plot(ts,x(:,3))
plot(ts,x(:,4))
legend('aj','bj','a2j','b2j')
xlabel('t')
ylabel('~Amplitude')

subplot(2,2,2),
plot(ts,H)
legend('Hamiltonian')
xlabel('t')
ylabel('~Energy')

subplot(2,2,3),
plot(x(:,1),x(:,2))
xlabel('aj')
ylabel('bj')

subplot(2,2,4),
plot(x(:,3),x(:,4))
xlabel('a2j')
ylabel('b2j')

%Define 1st order approximations

aj_first = @(t) (Aj*exp(1i*omega(j)*t)) + conj(Aj*exp(1i*omega(j)*t));
bj_first = @(t) (1i*g/omega(j))*Aj*exp(1i*omega(j)*t) + conj(1i*g/omega(j)*Aj*exp(1i*omega(j)*t));
a2j_first = @(t) my_alpha(2*j)*(Aj)^2*exp(2i*omega(j)*t) ...
    + conj(my_alpha(2*j)*(Aj)^2*exp(2i*omega(j)*t)) + my_gamma(2*j)*abs(Aj)^2;
b2j_first = @(t) 1i*my_beta(2*j)*(Aj)^2*exp(2i*omega(j)*t) + conj(1i*my_beta(2*j)*(Aj)^2*exp(2i*omega(j)*t));


%Plot numerical solution vs 1st order
figure()
subplot(2,2,1)
plot(ts,x(:,1))
hold on
plot(ts,aj_first(ts))
legend('numeric aj','1st order aj')

subplot(2,2,2)
plot(ts,x(:,2))
hold on
plot(ts,bj_first(ts))
legend('numeric bj','1st order bj')

subplot(2,2,3)
plot(ts,x(:,3))
hold on
plot(ts,a2j_first(ts))
legend('numeric a2j','2nd order a2j')

subplot(2,2,4)
plot(ts,x(:,4))
hold on
plot(ts,b2j_first(ts))
legend('numeric b2j','2nd order b2j')

%Functions

function wave_num = k(j)
wave_num = pi*j;
end

function tanh_kh = T(j)
h = 0.5; 
tanh_kh = tanh(abs(k(j))*h);
end

function freq = omega(j)
g=1;
freq = sqrt(g*k(j)*T(j));
end

function right_aj = RHS_aj(aj,bj,a2j,b2j,j)
g = 1;
linear = k(j)*T(j)*bj;
quad = 2*k(j)^2*(1-T(j)*T(2*j))*aj*b2j-k(j)^2*(1+T(j)^2)*a2j*bj;
cubic = -k(j)^3*T(j)*(3-2*T(j)*T(2*j))*aj^2*bj;
right_aj = linear+quad+cubic;
end

function right_bj = RHS_bj(aj,bj,a2j,b2j,j)
g = 1;
linear = -g*aj;
quad = -2*k(j)^2*(1-T(j)*T(2*j))*bj*b2j;
cubic = k(j)^3*T(j)*(3-2*T(j)*T(2*j))*aj*bj^2;
right_bj = linear+quad+cubic;
end

function right_a2j = RHS_a2j(aj,bj,a2j,b2j,j)
g = 1;
linear = k(2*j)*T(2*j)*b2j;
quad = 2*k(j)^2*(1-T(j)*T(2*j))*aj*bj;
right_a2j = linear+quad;
end

function right_b2j = RHS_b2j(aj,bj,a2j,b2j,j)
g = 1;
linear = -g*a2j;
quad = 1/2*k(j)^2*(1+T(j)^2)*bj^2;
right_b2j = linear+quad;
end

function slope = f(x,j)
aj = x(1,1);
bj = x(1,2);
a2j = x(1,3);
b2j = x(1,4);
slope = [RHS_aj(aj,bj,a2j,b2j,j) RHS_bj(aj,bj,a2j,b2j,j) RHS_a2j(aj,bj,a2j,b2j,j) RHS_b2j(aj,bj,a2j,b2j,j)];
end

function h3_up = h3(c1,c2,c3)
h3_up = -(k(c1)*k(c2) + abs(k(c1))*T(abs(c1))*abs(k(c2))*T(abs(c2)));
end

function H= hamiltonian(x,j)
g=1;
aj = x(1,1);
bj = x(1,2);
a2j = x(1,3);
b2j = x(1,4);
H = (g*aj^2+k(j)*T(j)*bj^2)+(g*a2j^2+k(2*j)*T(2*j)*b2j^2)-k(j)^2*(1+T(j)^2)*bj^2*a2j ...
    +4*k(j)^2*(1-T(j)*T(2*j))*b2j*bj*aj -k(j)^3*T(j)*(3-2*T(j)*T(2*j))*bj^2*aj^2;
end

function alpha = my_alpha(j)
alpha = k(j/2)*(3-T(j/2)^2)/(2*T(j/2)^3);
end

function beta = my_beta(j)
g=1;
beta = g*k(j/2)/omega(j/2)*(3+T(j/2))/(4*T(j/2)^3);
end

function gamma = my_gamma(j)
gamma = k(j/2)*(1+T(j)^2)/T(j);
end








