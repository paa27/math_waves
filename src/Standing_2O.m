%NL Standing wave (2nd-3rd order)
clear all
j = 1;
FFT = false;
wave_plot = true;
approx_plot = false;
plot_time = false;


dt=0.01; %Step size
N=3000; %# of total time steps
x=zeros(N,4);%initialze array
%x(1,:)=[0.0 sqrt(grav(2*j))*(1/k(j))^(3/2)*1/sqrt(2) 1/(2*k(j)) 0.0];
x(1,:)=[0.01 0. 0.000 0.0];
%Aj = x(1,1)/2 - 1i*x(1,2)*k(j)*T(j)/(2*omega(j,0));
%x(1,3) = my_alpha(2*j)*Aj^2+ conj(my_alpha(2*j)*Aj^2)...
%+ my_gamma(2*j)*abs(Aj)^2; %Initial conditions to set homogenous to 0
%x(1,4) = 1i*my_beta(2*j)*Aj^2 + conj(1i*my_beta(2*j)*Aj^2); %Initial conditions to set homogenous to 0
ts=0:dt:(N-1)*dt;
%xs=0:0.1:sqrt(2)*pi;
H=zeros(1,N);
H(1,1) = hamiltonian(x(1,:),j,0);

%standing = @(x,r) (2*cos(k(r).*x));

%Find 3rd order corrected resonance (deep water)

%res = @(x) (8*Aj^4.*x.^6+8.*Aj^2.*x.^4*(Aj^2-1)+x.^2.*(-8*Aj^2-2)+1);
%res_k = fzero(res,1/sqrt(2));

%RK4 method

for n=2:N
    ti = ts(n);
    g = grav(j,ti);
    x_old = x(n-1,:);
    
    k1 = dt.*f(x_old,j,ti);
    k2 = dt.*f(x_old+k1./2,j,ti+dt/2);
    k3 = dt.*f(x_old+k2./2,j,ti+dt/2);
    k4 = dt.*f(x_old+k3,j,ti+dt);
    
    x(n,:) = x_old + 1/6.*(k1+2.*k2+2.*k3+k4); 
    H(1,n) = hamiltonian(x(n,:),j,ti);
    
end

if plot_time == true
    
    [X,Y] = meshgrid(xs,ts); 
    standing1 = standing(xs,1);
    standing2 = standing(xs,2);
    [x1_up,x1_down] = envelope(x(:,1));
    [x3_up,x3_down] = envelope(x(:,3));
    my_sol1_down = envelope(standing1.*x(:,1));
    my_sol2_down = envelope(standing2.*x(:,3));
    my_sol1_up = standing1.*x1_up;
    my_sol2_up = standing2.*x3_up;
    figure()
    hold on
    for l=1:length(ts)/100
    plot(xs,my_sol1_up(1+100*(l-1),:)+real(Aj)*(l-1))    
    end
    hold off
    figure()
    surfplot1_up = surf(X,Y,my_sol1_up+my_sol2_up);
    hold on
    %surfplot1_down = surf(X,Y,my_sol1_down+my_sol2_down);
    xlabel('pos')
    ylabel('time')
    set(surfplot1_up,'LineStyle','none')
    %set(surfplot1_down,'LineStyle','none')
    title('j+2j mode')
    figure()
    surfplot2_up = surf(X,Y,my_sol1_up);
    hold on
    %surfplot2_down = surf(X,Y,my_sol1_down);
    xlabel('pos')
    ylabel('time')
    set(surfplot2_up,'LineStyle','none')
    %set(surfplot2_down,'LineStyle','none')
    title('j mode')
    figure()
    surfplot3_up = surf(X,Y,-my_sol2_up);
    hold on
    %surfplot3_down = surf(X,Y,my_sol2_down);
    xlabel('pos')
    ylabel('time')
    set(surfplot3_up,'LineStyle','none')
    %set(surfplot3_down,'LineStyle','none')
    title('2j mode')
    
figure()
plot(ts,x1_up,'r')
hold on
plot(ts,x1_down,'r')
plot(ts,x3_up,'b')
plot(ts,x3_down,'b')
legend('aj_{upper}','aj_{lower}','a2j_{upper}','a2j_{lower}')
xlabel('t')
ylabel('~Amplitude')
    
end



if wave_plot == true
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

fprintf('time_max2          2max        time_min1      1min         \n')

[a1_min, i_a1_min] = min(abs(x(:,1)));
[a2_max, i_a2_max] = max(abs(x(:,3)));
a1_min_t = ts(i_a1_min);
a2_max_t = ts(i_a2_max);

fprintf(' %f    	%f                   %f          %f  \n',a2_max_t,a2_max,a1_min_t,a1_min);

end

%Define 1st order approximations

if approx_plot == true

aj_first = @(t) (Aj*exp(1i*omega(j)*t)) + conj(Aj*exp(1i*omega(j)*t));
bj_first = @(t) (1i*g/omega(j))*Aj*exp(1i*omega(j)*t) + conj(1i*g/omega(j)*Aj*exp(1i*omega(j)*t));
a2j_first = @(t) my_alpha(2*j)*(Aj)^2*exp(2i*omega(j)*t) ...
    + conj(my_alpha(2*j)*(Aj)^2*exp(2i*omega(j)*t)) + my_gamma(2*j)*abs(Aj)^2;
b2j_first = @(t) 1i*my_beta(2*j)*(Aj)^2*exp(2i*omega(j)*t) ...
    + conj(1i*my_beta(2*j)*(Aj)^2*exp(2i*omega(j)*t));


%Plot numerical solution vs 1st order
figure()
subplot(2,2,1)
plot(ts,x(:,1))
hold on
plot(ts,aj_first(ts))
legend('numerical aj','1st order aj')

subplot(2,2,2)
plot(ts,x(:,2))
hold on
plot(ts,bj_first(ts))
legend('numerical bj','1st order bj')

subplot(2,2,3)
plot(ts,x(:,3))
hold on
plot(ts,a2j_first(ts))
legend('numerical a2j','2nd order a2j')

subplot(2,2,4)
plot(ts,x(:,4))
hold on
plot(ts,b2j_first(ts))
legend('numerical b2j','2nd order b2j')

end

%Check numerical frequency vs. non linear correction

if FFT == true

y = x(:,1);
my_start = 0;
my_end = 0;

for s=1:length(y)
    if y(s,1)<0 && y(s+1,1)>0
        my_start = s;
        break
    end
end

for s=my_start+10000:length(y)
    if y(s,1)<0 && y(s+1,1)>0
    my_end = s;
    break
    end
end

y = y(my_start:my_end,1);
% Sampling frequency
Fs = 1/dt;
% Calculate fft
ydft = fft(y);
% Only take one side of the Fourier transform
ydft = 2*ydft(1:ceil((length(y)+1)/2));
% Calculate the frequencies
freq = 0:Fs/length(y):Fs/2;
% Normalise according to length of signal
ydft = ydft/(2*length(freq));
[M,I] = max(ydft);
numerical_omega = freq(I)*2*pi
corrected_omega = omega(j)*(1+((9-12*T(j)^2-3*T(j)^4-2*T(j)^6)/(4*T(j)^4))*k(j)^2*Aj^2)
figure,
subplot(2,1,1), plot(ts(my_start:my_end),y), xlabel('Time [s]')
subplot(2,1,2), bar(2*pi*freq,abs(ydft)), set(gca,'yscale','log'), xlabel('omega')
end

%Functions

function wave_num = k(j)
res_k = pi;
%res_k = 0.701894158272629;
%L = (sqrt(2)*pi);
wave_num = res_k*j;
end

function tanh_kh = T(j)
h = 3000; 
tanh_kh = tanh(abs(k(j))*h);
end

function gravity = grav(j,ti)
g0 = 1;
sigma0 = 1;
freq = sqrt(g0*(1+sigma0/g0*k(1)^2)*k(1)*T(1));
gravity = g0*(1+sigma0/g0*k(j)^2)+0.3*cos(2*freq*ti);
end


function freq = omega(j,ti)
g=grav(j,ti);
freq = sqrt(g*k(j)*T(j));
end

function right_aj = RHS_aj(aj,bj,a2j,b2j,j,ti)
g = grav(j,ti);
linear = k(j)*T(j)*bj;
quad = 2*k(j)^2*(1-T(j)*T(2*j))*aj*b2j-k(j)^2*(1+T(j)^2)*a2j*bj;
cubic = -k(j)^3*T(j)*(3-2*T(j)*T(2*j))*aj^2*bj;
right_aj = linear+quad+cubic;
end

function right_bj = RHS_bj(aj,bj,a2j,b2j,j,ti)
g = grav(j,ti);
sigma0 = 1;
linear = -g*aj;
quad = -2*k(j)^2*(1-T(j)*T(2*j))*bj*b2j;
cubic = k(j)^3*T(j)*(3-2*T(j)*T(2*j))*aj*bj^2+3/2*sigma0*k(j)^4*aj^3;
right_bj = linear+quad+cubic;
end

function right_a2j = RHS_a2j(aj,bj,a2j,b2j,j,ti)
g = grav(2*j,ti);
linear = k(2*j)*T(2*j)*b2j;
quad = 2*k(j)^2*(1-T(j)*T(2*j))*aj*bj;
right_a2j = linear+quad;
end

function right_b2j = RHS_b2j(aj,bj,a2j,b2j,j,ti)
g = grav(2*j,ti);
linear = -g*a2j;
quad = 1/2*k(j)^2*(1+T(j)^2)*bj^2;
right_b2j = linear+quad;
end

function slope = f(x,j,ti)
aj = x(1,1);
bj = x(1,2);
a2j = x(1,3);
b2j = x(1,4);
slope = [RHS_aj(aj,bj,a2j,b2j,j,ti) RHS_bj(aj,bj,a2j,b2j,j,ti) ...
    RHS_a2j(aj,bj,a2j,b2j,j,ti) RHS_b2j(aj,bj,a2j,b2j,j,ti)];
end

function h3_up = h3(c1,c2,c3)
h3_up = -(k(c1)*k(c2) + abs(k(c1))*T(abs(c1))*abs(k(c2))*T(abs(c2)));
end

function H= hamiltonian(x,j,ti)
g1=grav(j,ti);
g2=grav(2*j,ti);
aj = x(1,1);
bj = x(1,2);
a2j = x(1,3);
b2j = x(1,4);

quad = (g1*aj^2+k(j)*T(j)*bj^2)+(g2*a2j^2+k(2*j)*T(2*j)*b2j^2);
cubic = -k(j)^2*(1+T(j)^2)*bj^2*a2j+4*k(j)^2*(1-T(j)*T(2*j))*b2j*bj*aj;
fourth = -k(j)^3*T(j)*(3-2*T(j)*T(2*j))*bj^2*aj^2;

H = quad + cubic;
end

function alpha = my_alpha(j)
alpha = k(j/2)*(3-T(j/2)^2)/(2*T(j/2)^3);
end

function beta = my_beta(j,ti)
g=grav(j,ti);
beta = g*k(j/2)/omega(j/2,ti)*(3+T(j/2))/(4*T(j/2)^3);
end

function gamma = my_gamma(j)
gamma = k(j/2)*(1+T(j)^2)/T(j);
end












