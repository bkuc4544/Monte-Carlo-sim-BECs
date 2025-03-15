%   li_sol - Program to solve the time dependent Schroedinger equation
%   using split operator technique


tic;

%  1 = w
%  1 = sqrt( hbar/(m*w) )
%  Using a trap frequency of about 50 hz this gives us a length scale
%  of about 5e-7 meters or 0.5 micron.
 
%length scale
L = 300; % units are quantum harmonic oscillator natural length
ngrid = 2^(14);
dz = L / ngrid;  % Coordinate step
dkz = 2*pi/L;    % Momentum step
kzmax = pi/dz;   % Maximum momentum
z = dz * (0:ngrid-1) - L/2; % Coordinate mesh
kz = dkz * (0:ngrid-1) - kzmax; % Momentum mesh
kzp = fftshift(kz); % New momentum mesh for FFT
T = 20; % max time, units of axial harmonic oscillator      20/ms
ntmax = 1000; % Number of time step
dt = T / ntmax; % Time step
dti = 1; % time step for imag time.
ntimax = 300; % max time step for imag time
 
%%%%% Hamiltonian %%%%%

Vharmonic = 0.5*(z).*(z); % Harmonic oscillator
DipBarrWidth = 100;  % Width of the green dipole barrier  
DipBarrAmp = 1.5*0.5*(DipBarrWidth^2);  % Amplitude of the green dipole barrier  (Currently set for second order cancelation)
VDipBarrier = DipBarrAmp*exp(-(z/DipBarrWidth).^2);
DipTrapWidth = 0.2;
DipTrapAmp = 0;
VDipTrap = -DipTrapAmp*exp(-(z/DipTrapWidth).^2);
RampSlope = 0;
Vramp = RampSlope*z;

V = Vharmonic + VDipBarrier ;
K = 0.5*kzp.^2 ; % Kinetic energy (momentum representation)
 
UV2 = exp(-1i * V * dt/2);  % Time evolution by potential energy(half time)
UK = exp(-1i * K * dt);  % Time evolution by kinetic energy

%50 micron width
 
%%%%% Initial wavefunction (ground state for SHO) %%%%%
sigma_z = 6;
z0 = 0;
psi_init = (exp(-((z-z0).^2)/(sigma_z^2)));
psi_initn = sum(psi_init.*conj(psi_init))*dz;
psi_init = psi_init/sqrt(psi_initn);
figure(1);
plot(z,psi_init.*(conj(psi_init)),'r');
hold on;

psi0i = psi_init;


plot(z,psi0i.*conj(psi0i),'b')
hold off;

%figure(2)
%imagesc(npsii)
 
 
% time evolution
ip = 0;
figure(3)
plot(V)
hold on


psi0 = psi0i;
CapDecayT = 0.5;  % On same scale as variable T
for ii=1:ntmax % for each time step
   if (mod(ii+9,10)== 0) 
       ip = ip+1;
       tray(ip) = ii*dt;
       npsi(ip,:) = psi0.*conj(psi0);
       psik = fft(psi0);
       apsik = psik.*conj(psik);
       napsik = sum(apsik)*dkz;
       KE(ip) = sum(K.*apsik/napsik)*dkz;
       
       zbar = sum((psi0.*z).*conj(psi0));
       varz = sum((psi0).*(z.*z).*conj(psi0));
       deltaz = varz^0.5;
       ExpectationValue(ip,:) = real(zbar);
       Uncertainty(ip,:) = real(deltaz);
   end  
   tt = ii*dt;
%    tray(ii)=tt;
   pnl = psi0.*conj(psi0);
   UV2_mf=exp(-1i * (V+VDipTrap*exp(-ii*dt/(CapDecayT))) * dt/2);
   %UV2_mf=exp(-i_imag * (Vd*(ii/ntmax)+Ng*psi0.*conj(psi0)) * dt/2);
   psi0 = UV2_mf .* ifft(UK .* fft(UV2_mf .* psi0)); % Split operator scheme

%    psi0nn(ii) = sum((psi0.*conj(psi0)))*dz;
%    psik = fft(psi0);
%    apsik = psik.*conj(psik);
%    napsik = sum(apsik)*dkz;
%    KE(ii) = sum(K.*apsik/napsik)*dkz;
%    apsi0 = psi0.*conj(psi0);
%    xbar(ii) = sum(z.*apsi0)*dz;
%    PE(ii) = sum(V.*apsi0*dz);
%    ME(ii) = sum(Ng*apsi0.*apsi0*dz/2);
%    EE(ii) = KE(ii)+PE(ii)+ME(ii);
end
 
%xbar = npsi*z'*dz/xbar0;
%PE = npsi*V'*dz;
%KE2 = apsik*K'/napsik*dkz;
%ME = Ng*sum(npsi.^2,2)*dz/2;
%EE = KE+PE'+ME';

%figure(7)
%plot(ExpectationValue)
%figure(8)
%plot(Uncertainty)

figure(4)
imagesc(z,tray,npsi,[0 .1])
colormap gray
figure(5)
imagesc((sign(log10(npsi)+7)+1).*(log10(npsi)+7))
%figure(6)
%plot(tray,KE,tray,PE,tray,ME,tray,EE)
%hold on
%plot(tray,xbar,'k-')
figure(1)
hold on
plot(z, npsi(end,:),'g')


%Code for computing the width of the wave function
zbar = sum((psi0.*z).*conj(psi0));
varz = sum((psi0).*(z.*z).*conj(psi0));
deltaz = varz^0.5;

toc

