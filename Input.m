%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAM's OUTPUT: Jacobian of Phase Velocity
%
% PROGRAMMERS:
% Prabir Das and Tarun Naskar
%
% Efficient analytical partial derivatives of modal phase velocity with respect to layer parameters
% P Das, T Naskar
% Geophysical Journal International 240 (3), 2091-2110
%
% Last revision date: 11/02/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Input.m is a MATLAB script to calculate the Jacobian of phase 
% velocity w.r.t S-wave velocity, P-wave velocity, Density, and thickness. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clearvars;
close all;
warning off
%% INPUT
% Model parameters
Model = xlsread('Soil_Profile.xlsx',1); 
Vs = Model(1,:);                % shear wave velocity
Vp_poi = Model(2,:);            % Vp_poi can be P-wave velocity or poisson's ratio
rho = Model(3,:);               % density
h = Model(4,:);                 % layer thickness
% Frequency Limits and Rayleigh Wave Mode
f_min = 0;                      % minimum frequency
f_max = 40;                     % maximum frequency
df = 0.5;                       % frequency resolution
dc = 1;                         % velocity resolution
mode = 1;                       % mode number
%% NAVIGATE TO THE FUNCTION CODES
addpath('sub_codes');
%% JACOBIAN OF PHASE VELOCITY w.r.t S-WAVE VELOCITY, P-WAVE VELOCITY, DENSITY, THICKNESS
[f,dVr_dVs,dVr_dVp,dVr_dVrho,dVr_dVh] = Jacobian(Vs, Vp_poi, rho, h, mode, f_min, f_max, df, dc) ;
%% PLOT THE JACOBIAN
figure;
plot(f,dVr_dVs,'LineWidth',1.2) % Ploting Jacobian of phase velocity w.r.to S-wave velocity 
set(gca,'TickDir', 'out','fontsize',16,'FontName','Times New Roman')
xlabel('Frequency (Hz)','FontSize',18,'FontWeight','bold','FontName','Times New Roman')
ylabel('Jacobian w. r. to Vs','FontSize',18,'FontWeight','bold','FontName','Times New Roman')
legend('Larer 1','Larer 2','Larer 3','Larer 4','Larer 5','Larer 6')

figure;
plot(f,dVr_dVp,'LineWidth',1.2) % Ploting Jacobian of phase velocity w.r.to P-wave velocity 
xlabel('Frequency(Hz)');
set(gca,'TickDir', 'out','fontsize',16,'FontName','Times New Roman')
xlabel('Frequency (Hz)','FontSize',18,'FontWeight','bold','FontName','Times New Roman')
ylabel('Jacobian w. r. to Vp','FontSize',18,'FontWeight','bold','FontName','Times New Roman')
legend('Larer 1','Larer 2','Larer 3','Larer 4','Larer 5','Larer 6')

figure;
plot(f,dVr_dVrho,'LineWidth',1.2) % Ploting Jacobian of phase velocity w.r.to density
xlabel('Frequency(Hz)');
set(gca,'TickDir', 'out','fontsize',16,'FontName','Times New Roman')
xlabel('Frequency (Hz)','FontSize',18,'FontWeight','bold','FontName','Times New Roman')
ylabel('Jacobian w. r. to Density','FontSize',18,'FontWeight','bold','FontName','Times New Roman')
legend('Larer 1','Larer 2','Larer 3','Larer 4','Larer 5','Larer 6')

figure;
plot(f,dVr_dVh,'LineWidth',1.2) % Ploting Jacobian of phase velocity w.r.to layer thickness
xlabel('Frequency(Hz)');
set(gca,'TickDir', 'out','fontsize',16,'FontName','Times New Roman')
xlabel('Frequency (Hz)','FontSize',18,'FontWeight','bold','FontName','Times New Roman')
ylabel('Jacobian w. r. to Layer-thickness','FontSize',18,'FontWeight','bold','FontName','Times New Roman')
legend('Larer 1','Larer 2','Larer 3','Larer 4','Larer 5','Larer 6')