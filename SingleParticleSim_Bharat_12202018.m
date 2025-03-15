clear all;
clc;
tic;

%Length of Simulation
T1 = 2*pi;  % total time in quantum harmonic oscillator time 1/w
T2 = 0.025*pi;
T3 = 20*pi;
T = T1+T2+T3;
nt = 2006;
dt = T/nt;
t_tot = dt*[0:nt-1];

%Problem is scaled to quantum harmonic oscillator length sqrt(hbar/mw)
%for w = 20hz*2*pi
% Mass = 100 Li7 atoms
L = 40;  
n = 2^12;
dx = L/n;
x = dx*[0:n-1]-L/2;

dk = 2*pi/L;
kmax = pi/dx;
k = dk*[0:n-1]-kmax;
k = fftshift(k);

% Dipole Barrier created by the green beam
DBoffset = 0;
DipBarrWidth = 135;  % Width of the green dipole barrier  
DipBarrAmp = 0.99 * 0.5 * (DipBarrWidth^2);  % Amplitude of the green dipole barrier  (Currently set for second order cancelation)
VDipBarrier =DipBarrAmp*exp(-((x-DBoffset)/DipBarrWidth).^2);

%Hamiltonian
KE = 0.5*k.^2;

% Unitary time evolution terms
UL = exp(-1i*KE*dt);

%Initialize the Wave Function
xo = 0;
s = 1;     % Equivalent to 0.5 micron width initially for the wavefunction
psi0 = exp(-(x-xo).^2/(2*s^2));
psi0 = psi0/sqrt(sum(psi0.*conj(psi0))*dx);  %normalization


%%Main Loop%%
img = zeros(nt,n);
psi = psi0;
for ii = 1:nt
    
    if (ii <= nt*(T1/T)) || (ii > nt*(T1+T2)/T)
        V = 0.5*x.^2 + VDipBarrier - DipBarrAmp;
    else
        V = 0.5*x.^2 + 0.9*(VDipBarrier - DipBarrAmp);
    end
    
    UN = exp(-1i*V*dt);
    img(ii,:) = psi.*conj(psi);
    psi = ifft( UL.*fft( UN.*psi ));  % Split Step
end

figure(1)
hold on;
plot(x,V)       % This will only plot V with the Dipole Barrier
xlabel('Position')
ylabel('Potential')
title('Potential flattening in the presence of a green dipole barrier')
hold off;

figure(2)
hold on;
plot(x,psi.*conj(psi),x,psi0.*conj(psi0))
legend('|\psi(4.025\pi)|^2', '|\psi(0)|^2')
xlabel('Position')
ylabel('|\psi|^2')
title('Probability')
hold off;

figure(3)
imagesc(x,t_tot,img)

toc

