% template for main.m,
% usage:copy to new dir, rename to main.m, change parameters
clear all;clc;close all;
%*******configuration**********
conf_file();
runtime=100e-12;
tstep=1e-15;
savestep=10;
enable_Hk_x=1;%1(0):enable(disable) hard axis anisotropy
%**********paramaters**********
%constant
constantfile();
%dimensions
LFL=50e-9;WFL=50e-9;tFL=0.6e-9;
LHM=LFL*1.1;WHM=WFL*1.1;tHM=2e-9;
%known parameters
Hk_Fe_x=1e-3/(gam*1e-12);%[T] w_parallel = 1e-3[THz]
Hk_Gd_x=Hk_Fe_x;%
if enable_Hk_x
    Hk_Fe_z=-0.023/(gam*1e-12);%[T] w_perp = 0.023[THz]
    Hk_Gd_z=Hk_Fe_z;%[T]
else
    Hk_Fe_z=0;%[T]
    Hk_Gd_z=0;%[T]
end
Hex_Fe=27.4/(gam*1e-12);%[T] w_E = 27.4[THz]
Hex_Gd=Hex_Fe;%[T]3meV->3meV/6.9838\mu_B=7.4204
alp=0.0068;
%unknown parameters
Ms_Fe=1149;%[emu/cm3]=1149e3 A/m, see THz osc paper ->d=0.4 nm, \mu_Fe=0.4^3*1149e-24=7.9293\mu_B
Ms_Gd=Ms_Fe;%[emu/cm3]=1012e3 A/m, see THz osc paper ->d=0.4 nm, \mu_Gd=0.4^3*1012e-24=6.9838\mu_B
Hext=[0,0,0];
%% STT parameters
jc_STT=-250e9;%[A/m2]
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
inita_theta=85/180*pi;inita_phi=0;%first sublattice in FL
initb_theta=inita_theta+pi;initb_phi=0;%second sublattice in FL
ma_init=[sin(inita_theta)*cos(inita_phi),sin(inita_theta)*sin(inita_phi),cos(inita_theta)];
mb_init=[sin(initb_theta)*cos(initb_phi),sin(initb_theta)*sin(initb_phi),cos(initb_theta)];
mmmPL=[1,0,0];
%% others
TT=300;%[K]
if dipolee
    %to do
else
    K12Dipole=zeros(3,3);
end
%Dx=0.01968237864387906;Dy=0.01968237864387906;
%Dz=0.960635227939411;%from online calculator
Dx=0;Dy=0;Dz=0;
Demag_=[Dx,0,0;0,Dy,0;0,0,Dz];
%% calc
totstep=round(runtime/tstep);
%**********dynamics**********
rk4_4llg();
%save('final.mat')
if (0)
    subplot(2,1,1);
    plot(tt*1e12,mmx(:,1),tt*1e12,mmy(:,1),tt*1e12,mmz(:,1),'linewidth',2);
    xlabel('time(ps)');ylabel('m');
    legend('mx','my','mz')
    subplot(2,1,2);
    plot(tt*1e12,mmx(:,2),tt*1e12,mmy(:,2),tt*1e12,mmz(:,2),'linewidth',2);
    xlabel('time(ps)');ylabel('m');
    legend('mx','my','mz')
end
if (1)
    figure
    plot(tt*1e12,mmx(:,1),tt*1e12,mmy(:,1),tt*1e12,mmz(:,1),'linewidth',2);
    xlabel('time(ps)');ylabel('m');
    legend('mx','my','mz')
end
if (0)
    load('relax_m1m2.mat')
    figure
    plot(tt*1e12,mmz(:,1)','-b','LineWidth',2);
    hold on
    clear all
    load('../AFMml_new/relax_ml.mat');
    plot(t*1e12,mmAz','-r','LineWidth',1);
    xlabel('time(ps)','fontsize',15);ylabel('my','fontsize',15)
    %xlim([0,15]);ylim([-1.05,1.05]);
    set(gca,'fontsize',20)
    legend('1fs','5fs')
end
