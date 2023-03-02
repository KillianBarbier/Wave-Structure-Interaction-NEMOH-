% TP WSINT - TEMER
% Template script 

% Alternatively, use run section and %%
% READ CAREFULLY THE SCRIPT AND CHECK THE values !!

%% Clear
clear variables
close all

%% INPUTS

% environment
depth =	250 ; %[m]  depth
nw =    100 ; % number of frequencies
w =     linspace(0.01,2, nw); % list of angular frequencies (rad/s)
dir =   0.0; % [deg] incident wave direction

% Target number of mesh panels
nb_panels = 150;

%======== END INPUTS===========%

%% PART 1. MESH

% Coarse description of the geometry
R = 7.2;        % radius of the cylinder at waterline
z_bottom = -78; % 
zG = -55;       % reference position for the computation of the Hydrodynamic database (here: position of the COG)


% Description of r, z and n for axiMesh_modified
%  Warning ! z(i) must be greater than z(i+1)
r = [R,R,0]        ; % array of radial coordinates 
z = [0,-78,-78]        ; % array of vertical coordinates

n=length(r);

% Target properties for the mesh
ntheta=ceil (sqrt(nb_panels)); % angular discretization. Different could be used !!

% Name of the save directory
dirname = "My_SPAR_"+num2str(nb_panels)+"_panels";

% Meshing function : axiMesh
[Mass,Inertia,KH,XB,YB,ZB] = axiMesh_modified(r,z,n,ntheta, zG, nb_panels, dirname);

%% PART 2. RUNNING NEMOH

% Call the Nemoh solver
[A,B,Fe]=Nemoh(w, dir, depth); %

% Save results into a Matlab .mat file in the save repertory
save(dirname+'/save.mat', 'A', 'B', 'w', 'Fe', 'KH', 'nb_panels');

%%  FIGURES : ADDED MASS, DAMPING, ...

% TODO: to complete

% Excitation force figures : modulus and phase

%ADDED MASS PLOT

plot_A({'My_SPAR_150_panels'},[1,1])
plot_A({'My_SPAR_150_panels'},[3,3])
plot_A({'My_SPAR_150_panels'},[5,5])

%DAMPING PLOT

plot_B({'My_SPAR_150_panels'},[1,1])
plot_B({'My_SPAR_150_panels'},[3,3])
plot_B({'My_SPAR_150_panels'},[5,5])

figure(2)
%
subplot(3,2,1)
plot(w, abs(Fe(:,1)))
ylabel('Fe_x')
title('Modulus')
%
subplot(3,2,2)
plot(w, phase(Fe(:,1)))
title('Phase')
%
subplot(3,2,3)
plot(w, abs(Fe(:,2)))
ylabel('Fe_y')
%
subplot(3,2,4)
plot(w, phase(Fe(:,2)))
%
subplot(3,2,5)
plot(w, abs(Fe(:,3)))
ylabel('Fe_z')
xlabel('Angular frequency (rad/s)')
%
subplot(3,2,6)
plot(w, phase(Fe(:,3)))
%
xlabel('Angular frequency (rad/s)')
sgtitle('Excitation force')



%% Mesh convergence

%figure
%load('My_SPAR_50_panels/save.mat')
%plot(...
%
nb_panels=[50,100,200,250];
for i=1:length(nb_panels)
       % Name of the save directory
       dirname = "My_SPAR_"+num2str(nb_panels(i))+"_panels";
       % Meshing function : axiMesh
       [Mass,Inertia,KH,XB,YB,ZB] = axiMesh_modified(r,z,n,ntheta, zG, nb_panels(i), dirname);
       % Call the Nemoh solver
       [A,B,Fe]=Nemoh(w, dir, depth); %

       % Save results into a Matlab .mat file in the save repertory
       save(dirname+'/save.mat', 'A', 'B', 'w', 'Fe', 'KH', 'nb_panels');
end

plot_Fe({'My_SPAR_50_panels','My_SPAR_100_panels','My_SPAR_150_panels','My_SPAR_200_panels','My_SPAR_250_panels',},[1,1])
plot_Fe({'My_SPAR_50_panels','My_SPAR_100_panels','My_SPAR_150_panels','My_SPAR_200_panels','My_SPAR_250_panels',},[3,3])
plot_Fe({'My_SPAR_50_panels','My_SPAR_100_panels','My_SPAR_150_panels','My_SPAR_200_panels','My_SPAR_250_panels',},[5,5])


%% OTHER METHOD 

% for i=1:length(nb_panels)
% % Target properties for the mesh
% ntheta=ceil (sqrt(nb_panels)); % angular discretization. Different could be used !!
% 
% % Name of the save directory
% dirname = "My_SPAR_"+num2str(nb_panels)+"_panels";
% 
% % Meshing function : axiMesh
% [Mass,Inertia,KH,XB,YB,ZB] = axiMesh_modified(r,z,n,ntheta, zG, nb_panels, dirname);
% % Call the Nemoh solver
% [A,B,Fe]=Nemoh(w, dir, depth); %
% 
% A1=A([1,3,5],[1,3,5],;)
% B1=B([1,3,5],[1,3,5],;)
% A2=vecnorm(A1);
% A3(:,i)=vecnorm(A2);
% B2=vecnorm(B1);
% B3(:,i)=vecnorm(B2);
% Fe1(:,i)=abs(Fe(:,1));
% Fe3(:,i)=abs(Fe(:,3));
% Fe5(:,i)=abs(Fe(:,5));
% Fe_overall=[Fe1,Fe3,Fe5];
% Fe7=vecnorm(Fe_overall);
% 
% % Save results into a Matlab .mat file in the save repertory
% save(dirname+'/save.mat', 'A', 'B', 'w', 'Fe', 'KH', 'nb_panels');
% 
% 
% end
% 
% save("A3.mat","A3")
% save("B3.mat","B3")
% save("Fe1.mat","Fe1")
% save("Fe3.mat","Fe3")
% save("Fe5.mat","Fe5")
% 
% load A3.mat
% load B3.mat
% load Fe1.mat
% load Fe3.mat
% load Fe5.mat

        


%% CALCULATE MOTIONS RAOs
addpath('plot')
nw=40
X_RAO = zeros(nw, 6);
X_RAO_moorings= zeros(nw, 6);
% TODO: IMPLEMENT RAO COMPUTATION HERE
i_xx=2.5*10^10;
i_yy=2.5*10^10;
i_zz=10^9;
km11=1.2*10^5;
km22=1.2*10^5;
%mass matrix
mass_matrix=[Mass 0 0 0 0 0;
               0 Mass 0 0 0 0;
               0 0 Mass 0 0 0;
               0 0 0 i_xx 0 0;
               0 0 0 0 i_yy 0;
               0 0 0 0 0 i_zz];
% Mooring stiffness matrix
KA=[km11 0 0 0 0 0;
    0 km22 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0];

for i=1:length(w)
    bdy_mtn=-(mass_matrix+A(:,:,i))*w(i)^2 + KH + KA - 1i*w(i)*B(:,:,i);
    X_RAO_moorings(i,:)=Fe(i,:)/bdy_mtn;
end
% PLOT MOTIONS RAOs with moorings

plot_RAO(w, X_RAO_moorings,'moorings');

%% PLOT without moorings

for i=1:length(w)
    bdy_mtn=-(mass_matrix+A(:,:,i))*w(i)^2 + KH  - 1i*w(i)*B(:,:,i);
    X_RAO(i,:)=Fe(i,:)/bdy_mtn;
end

plot_RAO(w, X_RAO);

%% Calculation eigen frequencies 

for i=1:6
    w_egn(i)=sqrt((KH(i,i)+KA(i,i))/(A(i,i)+mass_matrix(i,i)));
end

w_egn;
