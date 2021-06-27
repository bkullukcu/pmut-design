clc 
close all
syms X
syms x
%% Inputs
E3 = 2.450 * 10^9 ; %[Pa] Piezoelectric Young's Modulus 
v3 = 0.340 ;% Piezo Layer Poisson's Ratio
s11 = 1 / E3 ; % Compliance 1
s12 = -v3/E3 ; % Compliance 2
d31 = - 6 * 10^-12 ; % [C/N] 
e31f = d31/(s11+s12) ; % Effective Thin Film Piezoelectric Coeff. for PVDF
e31 = 0.024 ; % [C/m^2] Transverse Piezoelectric Coeff of PVDF
e33 = 9.5; % @ 233,6 kHz, 25 Celcius degrees , PVDF Relative Permitivity
e0 = 8.85*10^-12; % [F/m]
gama = 0.7; % Electrode radius to membrane radius ratio
a = 130*10^-6 ; % Plate radius [m]
area = pi*a^2 ; % Membrane Area
effarea = 1/3*area; %Effective Area of the Baffled Piston
ate = a*gama ; % Top Electrode radius [m]
n = 3 ; % Layer Number
phi = @(x) (1-x^2)^2 ; % Deflection shape function
phi1 =  diff((1-x^2)^2);  % Deflection shape function 1st Diff
phi2 = diff(diff((1-x^2)^2)) ; % Deflection shape function 2nd Diff
%modal = @(x) (((1-x^2)^2)^2)*x; % 
Z0 = 419; % [Rayl] Characteristic Impedance of Air
lambda00 = sqrt(10.22); % Eigenvalue of fundamental 01 vibration mode
v = 0.34 ; % Effective Poisson's Ratio for Composite Plate
aeff = a*((sqrt(3))/3) ; % [m] Effective Radius of the piston
t3 = 0.50*10^-6 ; % [m] Piezoelectric Layer Thickness
t2 = 1.50 * 10 ^-5 ;% [m] Structural Layer Thickness
t1 = 0.20*10^-6 ;%[m] Bottom electrode
v1 = 0.340 ;% Structural Layer Poisson's Ratio
v2 = 0.290;% Bottom electrode Poisson's Ratio
E1 = 2.500 * 10^9 ; %[Pa] Polyimide Young's Modulus
E2 = 27.8 * 10^9; %[Pa] Silver Young's modulus
z3 = t1+t2+(t3/2) ;%1.95 * 10^-6 ; %[m] Piezo layer central axis
z1 = t1/2 ;% 0.75 * 10^-6 ; %[m] Structural layer central axis
z2 = t1 + t2/2 ;%0.10 * 10^-6 ; %[m] Bottom electrode central axis
h3 = 2.20*10^-6 ; % [m] Top of Piezo layer
h1 = 1.50*10^-6 ; % [m] Top of Structural Layer
h2 = 1.70*10^-6 ; % [m] Top of Bottom electrode
ro3= 1780 ; % [kg/m^3] PVDF Density
ro1 = 1420 ; % [kg/m^3] PI Density
ro2 = 2700 ; % [kg/m^3] Aluminum Density
rote = 2700 ; % [kg/m^3] Aluminum Density 
roair = 1.225 ; % [kg/m^3] Air Density
tte = 0.2 * 10^-6 ; % [m] Top Electrode Thickness
vin = 10 ; % 10 [Volts] Input Voltage
r = 1*10^-2 ; % [m] Distance from the baffled piston
c = 343 ; % [m/s] Speed of sound in air medium
alpha = 0.9 ; % [Nepers/m] Absorbtion coefficient of air, %50 humidity, 20 Degrees Celcius,200 kHz
l = 250 * 10^-6 ; % [m] Tube length
plate1 = E1/(1-v1^2); % Plate Modulus for PI
plate2 = E2/(1-v2^2); % Plate Modulus for Bottom Electrode
plate3 = E3/(1-v3^2); % Plate Modulus for PVDF
w = 1*10^3 ; % [Hz] AC Frequency 
%% Equations
%% Electrical Domain
c0 = (e33*e0*pi*ate^2)/t3; % Capacitance
%% Mechanical Domain
zp = (((t1*z1*E1)+(t2*z2*E2)+(t3*z3*E3))/((t1*E1)+(t2*E2)+(t3*E3))); % Location of the neutral plane
zdist = z3 - zp; % Piezo Mid-Layer Distance to Neutral Axis
dist1 = h1 - zp ; % PI Distance to Neutral Axis
dist2 = h2 - zp ; % Bottom Electrode Distance to Neutral Axis
dist3 = h3 - zp ; % PVDF Distance to Neutral Axis
D = (1/3)*((plate1*(dist1^3))+(plate2*((dist2^3)-(dist1^3)))+(plate3*((dist3^3)-(dist2^3)))); % Flexural Rigidity
Ie = 5 ; % Strain Energy Integral
Im = -1.515 ; % Piezoelectric Coupling Integral
Km = ((a^2)/(2*pi*D*Ie))^-1; % 1/Mechanical Compliance
emcoupling = 2*Im*e31f*zdist; % Electromechanical Efficiency
w0 = (1/Km)*emcoupling*vin ; % Static Displacement
mdisk = ((ro1*t1+ro2*t2+ro3*t3)*pi*a^2)+rote*tte*pi*ate^2; % Total mass of the disk
Mm = mdisk* 0.2 ; % Modal Mass
%% Acoustics
%Ztube = (Z0)/(pi*a^2); % Impedance of acoustic transmisson line !!
mueff = (ro1*t1)+(ro2*t2)+(ro3*t3); % Effective Mass Per Unit Area
wn = ((lambda00/a)^2)*sqrt(D/mueff) ; % Natural Freq.
lambda = c/wn; % Wavelength
k = 2*pi/lambda ; % Wavenumber
Zac = ((Z0/(pi*aeff^2))/(1-(besselj(1,2*aeff)/aeff)+(1j*StruveH1(2*aeff)/aeff))); % Acoustic Impedance on front side of the transducer
Ztotal = real(Zac);
bair = ((pi*aeff^2)^2)*ZtotalRe; % Total air damping
bw = (bair)/(2*pi*Mm); % Bandwidth
bwp = (1-((wn-bw)/wn))*100; % Percent Bandwidth
Q = (sqrt(Km*Mm))/bair; % Quality Factor
%{
%Zth = ((Z0/(pi*a^2))/(1-(besselj(1,2*a)/a)+(1j*StruveH1(2*a)/a))); % Acoustic Impedance on open end of the tube !!
Rth = (Zth - Ztube) / (Zth + Ztube); % Reflection Coeff.
%Zback = Ztube*((exp(1i*k*l)+Rth*exp(-1i*k*l))/(exp(1i*k*l)-Rth*exp(-1i*k*l))) ; % Impedance of the back side of the plate !!
Tth = 2*Zth/Zth+Ztube ; % Transmisson Coefficients
Ra = real(Zac)*(effarea/emcoupling)^2; % Real Part of Acoustic Impedance
%Ztotal = Zfront ;%+ Zback ; 
%ZtotalRe = real(Ztotal);
Pin = (emcoupling/(pi*aeff^2)); % Input Pressure
zfm = Zac+(1/(pi*aeff^2)^2)*(1i*2*pi*w0*Mm+(Km/(1i*2*pi*w0))) ; 
Tp = (2*zfm)/(zfm + Ztube); 
Peq = (Pin/cos(k*l))*(-(1i*Ztube*cot(k*l))/(zfm-1i*Ztube*cot(k*l))) ; % Thevenin Eqv. Pressure
Gtx = Peq/Pin ;
Rp = (zfm - Ztube)/(zfm + Ztube);
Zeq = Ztube*((exp(1i*k*l)+Rp*exp(-1i*k*l))/(exp(1i*k*l)-Rp*exp(-1i*k*l)));
vel = abs((emcoupling/(pi^2)*(a^2)*(aeff^2))*(Gtx/Zth+Zeq)) ; % Particle Velocity
Pdist = ((pi*(a^2)*Z0*vel)/(lambda*r))*exp(-alpha*r); % Far field on axis pressure of a baffled piston at distance r
Grx = (-1i*Ztube*cot(k*l))/((Zth-1i*Ztube*cot(k*l))*cos(k*l));%e31/(e33*e0); % Receive Sensitivity
Zc0 = -1i*(1/w*c0); % Impedance of the capacitance
Srxb = (pi*aeff^2/emcoupling)*Grx*(Zc0/(zfm+Zback+Zc0));
Srxf = (pi*aeff^2/emcoupling)*(Zc0/(zfm+Zback+Zc0));
Preceive = Peq/Grx ; % Receive Pressure
%}