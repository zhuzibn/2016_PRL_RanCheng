% reproduce fig.1 of 2016-Terahertz Antiferromagnetic Spin Hall Nano-Oscillator-PRL-Ran Cheng
clear all;clc;close all
%% parameters obtained from ref.30, in the unit of 2*pi*THz
wperp=0.023;
wpara=0.001;
alph=0.0068;
wE=27.4;

ws=linspace(0,1,100)*wperp;
wp=wE*(1i*alph+sqrt((wperp+2*wpara+sqrt(wperp^2-4*ws.^2))./wE-alph^2));%optical mode
wn=wE*(1i*alph+sqrt((wperp+2*wpara-sqrt(wperp^2-4*ws.^2))./wE-alph^2));%acoustic mode
figure;hold on
plot(ws/wperp,real(wp),'-r')
plot(ws/wperp,real(wn),'-b')
plot(ws/wperp,imag(wp),'--k')
plot(ws/wperp,imag(wn),'--k')
xlim([0 1]);%ylim([0 1.2])
xlabel('{\omega}_s/{\omega}_{\perp}');ylabel('{\omega}(2{\pi}THz)')
%wsth=sqrt(wperp^2/4+alph^2*(2*wpara+wperp)*wE); %Eq.3
%save('fig1.mat')