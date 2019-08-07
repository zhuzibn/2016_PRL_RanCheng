% numerical solve Eq. S4 of 2016-Terahertz Antiferromagnetic Spin Hall Nano-Oscillator-PRL-Ran Cheng
clear all;clc;close all
%% parameters obtained from ref.30, in the unit of 2*pi*THz
approx=2;
%0: no approximation 1: approximation 2: complete expression, refer to
%2016-PRL-Ran Cheng.pdf for definition 
wperp=0.023;
wpara=0.001;
alph=0.0068;
wE=27.4;
ws_=linspace(0,1,100)*wperp;
szws=size(ws_,2);
w_=zeros(szws,4);
syms w
for ct=1:szws
    ws=ws_(ct);
    switch approx
        case 0
            tmp=vpasolve((2*wE*ws+wpara*ws)*(-2*wE*ws+(wpara+wperp)*ws)-(2*wE*(wpara+wperp+1i*alph*w)-w^2+wpara*(wpara+wperp+1i*alph*w))*(2*wE*(wpara+1i*alph*w)-w^2+(wpara+wperp)*(wpara+1i*alph*w))==0,w);
        case 1
            tmp=vpasolve((w^2-2*wE*(wpara+wperp+1i*alph*w))*(w^2-2*wE*(wpara+1i*alph*w))+(2*wE*ws)^2==0,w);
        case 2
            tmp=vpasolve((-2*(wE+wpara+1i*alph*w)*ws)*(2*(wE+wpara+wperp+1i*alph*w)*ws)-...
                (w^2+ws^2-(wpara+1i*alph*w)*(2*wE+wpara+wperp+1i*alph*w))*(w^2+ws^2-(wpara+wperp+1i*alph*w)*(2*wE+wpara+1i*alph*w))==0,w);
    end
    w_(ct,1)=tmp(1);
    w_(ct,2)=tmp(2);
    w_(ct,3)=tmp(3);
    w_(ct,4)=tmp(4);
end
%save('EqS4.mat');
figure;hold on
plot(ws_/wperp,real(w_(:,3)),'-r')
plot(ws_/wperp,real(w_(:,4)),'-b')
plot(ws_/wperp,imag(w_(:,3)),'k*')
plot(ws_/wperp,imag(w_(:,4)),'k*')
