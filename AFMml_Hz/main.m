% [1]. 2016-Terahertz Antiferromagnetic Spin Hall Nano-Oscillator-PRL-Ran Cheng
clear all;clc;close all;
%*******configuration**********
conf_file();
runtime=40e-12;
tstep=1e-15;
savestep=10;
%**********paramaters**********
%constant
constantfile();
%dimensions
LFL=50e-9;WFL=50e-9;tFL=0.6e-9;
LHM=LFL*1.1;WHM=WFL*1.1;tHM=2e-9;
%known parameters
alp=0.0068;
w_para=0.001*1e12;%[Hz]
w_perp=0.023*1e12;%[Hz]
wE=27.4*1e12;%[Hz]
alp_NL=0.015;%feedback coefficient
%unknown parameters
Ms=1149;%[emu/cm3]=1e6 A/m
Hk=0;
Hext=[0,0,0];
%% STT parameters
jc_STT=-200e9;%[A/m2]
PolFL=0.4;%polarization of FL layer
PolSTT=[1,0,0];
if STT_FLT
facFLT_STT=0.2;%ratio of FLT over DLT
else
facFLT_STT=0;    
end
%% SOT parameters
thetaSH=0.2;
lambdaSF=5e-9;%spin diffusion length
polSOT=[0,1,0];%spin flux polarization
jc_SOT=0e10;%[A/m2]
if SOT_FLT
    facFLT_SHE=2;%ratio of FLT/DLT
else
    facFLT_SHE=0;
end
%% initial condition
init_theta1=85/180*pi;init_phi1=0;% first sublattice
init_theta2=init_theta1+180/180*pi;init_phi2=0;% second sublattice
m_init1=[sin(init_theta1)*cos(init_phi1),sin(init_theta1)*sin(init_phi1),cos(init_theta1)];
m_init2=[sin(init_theta2)*cos(init_phi2),sin(init_theta2)*sin(init_phi2),cos(init_theta2)];
mmmPL=[1,0,0];
%% others
TT=300;%[K]

if dipolee
    %to do
else
   K12Dipole=zeros(3,3); 
end

Dx=0;Dy=0;Dz=0;
Demag_=[Dx,0,0;0,Dy,0;0,0,Dz];
%% calc
totstep=round(runtime/tstep);
%**********dynamics**********
rk4_4llg();
%save('final.mat')
if(0)
figure;
plot(t*1e12,mmy,'linewidth',2)
xlabel('time(ps)');ylabel('m')   
figure;
plot(t*1e12,Lx,'linewidth',2)
xlabel('time(ps)');ylabel('L')  
end
if(0)
subplot(2,1,1);
plot(t*1e12,mmx,t*1e12,mmy,t*1e12,mmz,'linewidth',2)
xlabel('time(ps)');ylabel('m')
legend('mx','my','mz')
subplot(2,1,2);
plot(t*1e12,Lx,t*1e12,Ly,t*1e12,Lz,'linewidth',2)
xlabel('time(ps)');ylabel('L')
legend('Lx','Ly','Lz')
end
if(1)
figure;
plot(t*1e12,mmAx,t*1e12,mmAy,t*1e12,mmAz,'linewidth',2)
xlabel('time(ps)');ylabel('m_A')  
legend('mx','my','mz')
end
if (0)%compare ml_T & ml_Hz relaxation
    clear all;close all
    load('relax_ml_Hz.mat')
    figure
    plot(t*1e12,mmAx'','-b','LineWidth',2);
    hold on
    clear all
    load('../AFMml_T/relax_ml_T.mat');
    plot(t*1e12,mmAx','-r','LineWidth',1);
    xlabel('time(ps)','fontsize',15);ylabel('mz','fontsize',15)
    %xlim([0,15]);ylim([-1.05,1.05]);
    set(gca,'fontsize',20)
    legend('1fs','5fs')
end
if (0)%compare ml_T & ml_Hz Jc=-200e9
    clear all;close all
    load('ml_Hz_-200.mat')
    figure
    plot(t*1e12,mmAz'','-b','LineWidth',2);
    hold on
    clear all
    load('../AFMml_T/ml_T_-200.mat');
    plot(t*1e12,mmAz','-r','LineWidth',1);
    xlabel('time(ps)','fontsize',15);ylabel('mz','fontsize',15)
    %xlim([0,15]);ylim([-1.05,1.05]);
    set(gca,'fontsize',20)
    legend('1fs','5fs')
end







