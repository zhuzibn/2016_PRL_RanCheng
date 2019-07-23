clear all;close all;clc
addpath('D:\tmp\2016PRLRanCheng\AFMm1m2_LL')
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
if (0)% one level
    datrange=[-190,-200,-220:-20:-500,-1000,-2000];
    %datrange=[-10,-20,-50,-100,-200,-500,-1000,-2000,-5000];
    szdate=size(datrange,2);
    mxFe_=zeros(szdate,1);mxGd_=zeros(szdate,1);
    angleFeGd_=zeros(szdate,1);
    Freqq_=zeros(szdate,1);
    angleFex_=zeros(szdate,1);
    for ctdat=1:szdate
        datname=sprintf('final%d.mat',datrange(ctdat));
        datname
        load(datname);
        if (0)
        subplot(2,1,1);
        plot(tt*1e12,mmx(:,1),tt*1e12,mmy(:,1),tt*1e12,mmz(:,1),'linewidth',2);
        legend('mx','my','mz')
        xlabel('time(ps)');ylabel('m');
        subplot(2,1,2);
        plot(tt*1e12,mmx(:,2),tt*1e12,mmy(:,2),tt*1e12,mmz(:,2),'linewidth',2);
        legend('mx','my','mz')
        xlabel('time(ps)');ylabel('m');
        end
        if (0)
        subplot(2,1,1);
        plot(tt*1e12,mmx(:,1),'linewidth',2);
        xlabel('time(ps)');ylabel('mx');
        subplot(2,1,2);
        plot(tt*1e12,mmx(:,2),'linewidth',2);
        xlabel('time(ps)');ylabel('mx');            
        end
        if (1)%FFT
            plotfft=1;
            nT0=size(mmy,1);
            runTime0=runtime*1e9;%[ns]
            y=mmy(:,1)';
            rminit=0.4;
            rmlast=0.01;
            Freqq_(ctdat)=FFT_module(nT0,runTime0,y,rminit,rmlast,plotfft);
        end
        close all;
        mxFe_(ctdat)=mmx(end-1,1);
        mxGd_(ctdat)=mmx(end-1,2);
        u=[mmx(end-1,1),mmy(end-1,1),mmz(end-1,1)];
        v=[mmx(end-1,2),mmy(end-1,2),mmz(end-1,2)];
        CosTheta = dot(u,v)/(norm(u)*norm(v));
        angleFeGd_(ctdat) = acosd(CosTheta);%angle bw two sublattices
        angleFex_(ctdat)=acosd(mxFe_(ctdat));
    end
end

if (0)%compare with previous result
    dat_=[1000,4000];
    datname1=sprintf('final%d',dat_(1));
    figure;hold on
    load(datname1);
    plot(tt*1e12,mmy(:,1),'-b','LineWidth',2);
    datname2=sprintf('final%d',dat_(2));
    load(datname2);
    plot(tt*1e12,mmy(:,1),'-r','LineWidth',1);
    xlabel('time(ps)','fontsize',15);ylabel('my','fontsize',15)
    %xlim([0,15]);ylim([-1.05,1.05]);
    set(gca,'fontsize',20)
    legend(datname1,datname2)
end
if (1)
    hbar=6.58211951440e-16;%eV.s
    gam=1.760859644e11;%rad/(s.T)
    tz=0.6e-9;%[m]
    msTM=1149e3;%[A/m]
    etaSTT=0.8;%spin hall angle
    JSTT_=[200:20:500]*1e9;
    BDSTTTM=hbar/2*etaSTT*JSTT_/(msTM*tz);
    wSTT=gam*BDSTTTM*1e-12;%[THz]
    Freqnum=[313.559322;347.4576271;381.3559322;406.779661;440.6779661;...
        474.5762712;500;533.8983051;567.7966102;601.6949153;627.1186441;...
        661.0169492;694.9152542;720.3389831;754.2372881;788.1355932]*1e-3;%THz
    wperp=0.023;
    wpara=0.001;%[THz];
    alph=0.0068;
    wE=27.4;%[THz]
    ws=linspace(0,1,100)*wperp;
    wn=wE*(1i*alph+sqrt((wperp+2*wpara-sqrt(wperp^2-4*ws.^2))./wE-alph^2));%acoustic mode
    figure;hold on
    plot(wSTT/wperp,Freqnum,'-*b')
    %plot(ws/wperp,real(wn)/2.4,'-r')
    plot(ws/wperp,real(wn),'-r')
    xlim([0 1]);%ylim([0 1.2])
    legend('numerical','analytical')
    xlabel('{\omega}_s/{\omega}_{\perp}');ylabel('{\omega}(2{\pi}THz)')
end