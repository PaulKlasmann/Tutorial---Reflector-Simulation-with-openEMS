% Parabolic Reflector Example
%
% A simple Ku band feed is used to illuminate a 0.34m parabolic reflector
% 
% Tested in Windows 7 with
%  - Octave 4.4.1
%  - openEMS v0.0.35
% Requires around 5 Gb of memory in addition to OS memory overheads
% Author: Paul Klasmann with a lot of help from Thorsten Liebig
% 01/12/2018
%
% In this script, the waveguide is the "horn" even though it is not flared.

close all
clear
clc
% Simulation constants
physical_constants; % calls a function that initializes c, ep0, mu0 and eta0
unit = 1e-3; % all lengths in mm

wg.radius  = 19/2;      % Horn radius
wg.length = 100;        % Horn length in z-direction
wg.thickness = 2;       % Horn's wall thickness
FL=155;                 % Reflectors focal length
R = 170;                % Reflector radius in mm
A = pi*(R*unit)^2;      % Aperture Area

% frequency range of interest Ku Satcom Band
f_start =  10e9;        % Start Frequency
f_stop  =  15e9;        % Stop Frequency
f0 = 12e9;              % Frequency of interest

x=[0:0.5:R];            % Vector to store x values for parabola

z=(((x).^2)/(4*FL))-FL; % Parabola curve definition in z axix and offset by FL
curve = [x;z];          % Store x and z coordinates in "curve"
depth = ((2*R)^2)/(16*FL);    % Calculate depth of reflector, depth=D^2/16FL
extrapoints = [R 0 0; -FL-2 -FL-2 0]; % Extra coords to make solid reflector
coords = [curve,extrapoints];         % Full coords for RotPoly stored here

plot(x,z);                            % Plot parabola to check
axis equal                            % Scale axis equally for aspect ratio 1:1
set(gca, "linewidth",2, "fontsize", 14 )
xlabel( 'Radius in x Direction (mm)', 'FontSize', 14 );
ylabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
title( 'Parabola Profile', 'FontSize', 16 );

% Size of the simulation box, should be at least lambda/4 at lowest frequency
SimBox = [380 380 500];               % Size in x, y and z directions

% Initialise the FDTD structure
FDTD = InitFDTD( 'NrTS', 5000, 'EndCriteria', 0.5e-3 ); % End criteria -33 dB
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % PML boundary
%BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % Simple MUR boundary
FDTD = SetBoundaryCond( FDTD, BC );   % Store FDTD and BC in FDTD structure

% Initialize the mesh with the "air-box" dimensions and key mesh lines for ports
max_res = c0 / (f_stop) / unit / 20; % cell size: lambda/20

mesh.x = [-SimBox(1)/2 -wg.radius 0 wg.radius SimBox(1)/2];
% create a smooth mesh between specified fixed mesh lines
mesh.x = SmoothMeshLines( mesh.x, max_res, 1.4);
mesh.y = mesh.x;                      % Structure is axis symmetric in x & y 
mesh.z = [-FL-50 -FL-2 -FL -FL+depth 0 wg.length-5 wg.length wg.length+1 wg.length+20];
mesh.z = SmoothMeshLines(mesh.z, max_res, 1.4);

%% Setup CSXCAD geometry & mesh
CSX = InitCSX();
CSX = DefineRectGrid( CSX, unit, mesh );
CSX = AddMetal(CSX, 'Parabola');
% x is normal direction, "coords" is array to plot shape in 2D, z is the rotational
% axis direction.
CSX = AddRotPoly(CSX,'Parabola',10,'x',coords,'z'); 

% Waveguide horn is a rotational polygon
% The open aperture of the horn is at z=0 (at the reflector's focal point)
CSX = AddMetal(CSX, 'Circular_Waveguide');
p(1,1) = wg.radius;   
p(2,1) = 0;
p(1,2) = wg.radius+wg.thickness;
p(2,2) = 0;
p(1,3) = wg.radius+wg.thickness;
p(2,3) = wg.length;
p(1,4) = wg.radius;
p(2,4) = wg.length;
CSX = AddRotPoly(CSX,'Circular_Waveguide',10,'x',p,'z');

% End cap to prevent the radiation coming out of the back of the "horn"
CSX = AddMetal(CSX,'Cap');
CSX = AddCylinder(CSX,'Cap',10,[0 0 wg.length],[0 0 wg.length+1],wg.radius+wg.thickness);
% End of model geometry

%% Apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start=[-wg.radius -wg.radius wg.length ];
stop =[+wg.radius +wg.radius wg.length-5 ];
[CSX, port] = AddCircWaveGuidePort( CSX, 0, 1, start, stop, wg.radius*unit, 'TE11', 0, 1);

% Dump box for Electric field at Phi=0 (vertical cut)
CSX = AddDump(CSX,'Et_V_dump', 'SubSampling', '4,4,4');
start=[0 -380/2 -FL];
stop =[0 380/2 wg.length+20];
CSX = AddBox(CSX,'Et_V_dump',0,start,stop);

% Dump box for Electric field at Phi=90 (horizontal cut)
CSX = AddDump(CSX,'Et_H_dump', 'SubSampling', '4,4,4');
start=[-380/2 0 -FL];
stop =[380/2 0 wg.length+20];
CSX = AddBox(CSX,'Et_H_dump',0,start,stop);

% nf2ff calc
start = [mesh.x(10) mesh.y(10) mesh.z(10)];
stop  = [mesh.x(end-11) mesh.y(end-11) mesh.z(end-11)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'Directions', [1 1 1 1 1 1], 'OptResolution', max_res*2);

% Prepare simulation folder
Sim_Path = 'tmp';
Sim_CSX = 'Parabola.xml';
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder
% Write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );
% Show the structure
%CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
CSXGeomPlot( [Sim_Path '/' Sim_CSX], ['--export-polydata-vtk=tmp'] );

% Run openEMS
%openEMS_opts = '--debug-PEC --no-simulation'; % uncomment to visualise mesh
RunOpenEMS( Sim_Path, Sim_CSX, '--numThreads=4');

% Postprocessing & do the plots
freq = linspace(f_start,f_stop,201);
port = calcPort(port, Sim_Path, freq);

Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;

% Plot reflection coefficient S11
figure
plot( freq/1e9, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
ylim([-40 0]);
set(gca, "linewidth",2, "fontsize", 14 )
grid on
title( 'Reflection Coefficient S_{11}', 'FontSize', 16 );
xlabel( 'Frequency (GHz)','FontSize', 14 );
ylabel( 'Reflection Coefficient |S_{11}| (dB)','FontSize', 14 );
drawnow

% NFFF plots

% calculate the far field at phi=0, 45 and at phi=90 degrees
thetaRange = (0:0.2:359) - 180;
disp( 'calculating far field at phi=[0 45 90] deg...' );
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, thetaRange*pi/180, [0 45 90]*pi/180);

Dlog=10*log10(nf2ff.Dmax);       % Calculate maximum Directivity in dB
G_a = 4*pi*A/(c0/f0)^2;          % Calculate theoretical gain for given aperture
e_a = nf2ff.Dmax/G_a;            % Calculate Efficiency

% Display some antenna parameters from above calculations
disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(Dlog) ' dBi'] );
disp( ['aperture efficiency: e_a = ' num2str(e_a*100) '%'] );

% Directivity
figure
plotFFdB(nf2ff,'xaxis','theta','param',[1 3]);
ylim([-30 40]);
xlim([-180 180]);
grid on
set(gca,"linewidth",2, "fontsize", 14, "XTick", -180:30:180, "YTick", -30:5:40 )
title( 'Farfield Directivity @ 12GHz','FontSize', 16 );
xlabel( 'Frequency (GHz)','FontSize', 14 );
ylabel( 'Directivity (dBi)','FontSize', 14 );
drawnow

% Plot Ludwig3 cross polar
plotFFcocx(nf2ff,'xaxis','theta','param',[2]);
ylim([-30 40]);
xlim([-180 180]);
grid on
set(gca,"linewidth",2, "fontsize", 14, "XTick", -180:30:180, "YTick", -30:5:40 )
drawnow

% Polar plot
figure
leg=[]; %legend
polarFF(nf2ff,'xaxis','theta','param',[1 3],'logscale',[-30 35], 'xtics', 12);
title( 'Farfield Directivity @ 12GHz','FontSize', 16 );
xlabel( 'Frequency (GHz)','FontSize', 14 );
ylabel( 'Directivity (dBi)','FontSize', 14 );
drawnow

%% Calculate 3D pattern
phiRange = sort( unique( [-180:5:-100 -100:2.5:-50 -50:1:50 50:2.5:100 100:5:180] ) );
thetaRange = sort( unique([ 0:1:50 50:2.:100 100:5:180 ]));

disp( 'calculating 3D far field...' );
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, thetaRange*pi/180, phiRange*pi/180, 'Verbose',2,'Outfile','nf2ff_3D.h5');

figure
plotFF3D( nf2ff, 'logscale', -40);        % plot 3D far field in dB

% Save far field in VTK to plot in ParaView
E_far_normalized = nf2ff.E_norm{1}/max(nf2ff.E_norm{1}(:));
DumpFF2VTK([Sim_Path '/Farfield.vtk'],E_far_normalized,thetaRange,phiRange,'scale', 0.005, 'logscale', -20, 'maxgain', 31);

