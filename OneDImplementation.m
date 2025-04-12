clear;clc;close all;
%Constants
ur = 2; er = 6;
c0 = 299792458; gigahertz = 1e9;

%Wavelength Grid Resolution
fmax = 1*gigahertz;
Nlamdda = 20;
nmax = sqrt(ur*er);
lambdamin = c0/fmax/nmax;
dlambda = lambdamin/Nlamdda;

%Structure Grid Resolution
d=30.48e-2;
Nd = 4;
deld = d/Nd;
    
%Grid res
dz = min(dlambda,deld);
Nc = ceil(d/dz);
cell = 10;%Spacer region
N = Nc+ 2*cell + 3;

nz1 = 13; %Initial position ---Change
nz2 = 83; %Final Position

%Material applied to these regions
UR(nz1:nz2) = ur;
ER(nz1:nz2) = er;


%Time step 
nbc = 1; %Reflection 
dt = nbc*dz/(2*c0);

%Source parameters
B = gigahertz; %Bandwith
tau = 1/(2*B);
t0 = 6*tau;
tprop = nmax*N*dz/c0; % time steps
T = 12*tau + 5*tprop;
STEPS = ceil(T/dt);

%Source Functions
nsrc = 1; ersrc = 1; ursrc = 1;
t = (0:STEPS-1)*dt;
delt = nsrc*dz/(2*c0) + dt/2;
A = -sqrt(ersrc/ursrc); %permittivity and permeability before material --->Change 
Esrc = exp(-((t-t0)/tau).^2);
Hsrc = A*exp(-((t-t0+delt)/tau).^2);

% Fourier Transforms
NFREQ = 100;
FREQ = linspace(0,1*gigahertz,NFREQ);
K = exp(-1i*2*pi*dt*FREQ);
REF = zeros(1,NFREQ);
TRN = zeros(1,NFREQ);
SRC = zeros(1,NFREQ);

Nz=N;
Ey = zeros(1, Nz);
Hx = zeros(1, Nz-1);
eps0 = 8.854187817e-12;
mu0 = 4*pi*1e-7;

% Fill in ER and UR arrays with background values if not already done
if length(ER) < N
    ER = [ones(1, N)];
    ER(nz1:nz2) = er;
end

if length(UR) < N
    UR = [ones(1, N)];
    UR(nz1:nz2) = ur;
end

% Define mEy and mHx
mEy = dt ./ (eps0 * ER);      % Length Nz
mHx = dt ./ (mu0 * UR(1:end-1)); % Length Nz-1


z = (0:Nz-1) * dz;  % Physical positions in meters

% Plot setup
figure;
hold on;
ePlot = plot(z, Ey, 'b', 'DisplayName', 'Ey');
mPlot = plot(z, [Hx, 0], 'r', 'DisplayName', 'Hx');  % Hx has Nz-1 points
legend();
ylim([-20, 20]);
xlabel('z-position (m)');
ylabel('Field Amplitude');
title('1D FDTD: Electric and Magnetic Fields')


h2 = 0; h1 = h2; e1 = h1; e2 = e1;

for T = 1:STEPS
    % Record H at boundary 
    h2 = h1;    h1 = Hx(1);

    % Update H
    for nz = 1:Nz-1
        Hx(nz) = Hx(nz) + mHx(nz)/dz * (Ey(nz+1) - Ey(nz));
    end
    Hx(end) = Hx(end) + mHx(end)/dz * (e2-Ey(end)); % Update H at Nz

    % Correct H at source (add source term)
    Hx(nz1) = Hx(nz1) + Hsrc(T);  % Apply the magnetic source

    % Record E at boundary
    e2 = e1;  e1 = Ey(Nz);

    % Update E
    for nz = 2:Nz-1
        Ey(nz) = Ey(nz) + mEy(nz)/dz * (Hx(nz) - Hx(nz-1));
    end
    Ey(1) = Ey(1) + mEy(1)/dz * (Hx(1) - h2); % Update E at 1

    % Correct E at source (add source term)
    Ey(nz1) = Ey(nz1) + Esrc(T);  % Apply the electric source

    % Update plots
    set(ePlot, 'YData', Ey);
    set(mPlot, 'YData', [Hx, 0]);  % Pad to match z axis length
    title(['1D FDTD: Ey and Hx | Time step: ', num2str(T), ' / ', num2str(STEPS)]);
    pause(0.0001);
end
