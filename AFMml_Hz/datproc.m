clear all;close all;clc
addpath('D:\tmp\2016PRLRanCheng\AFMml_Hz')
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
    datrange=[-100:-20:-500];
    %datrange=[-10,-20,-50,-100,-200,-500,-1000,-2000,-5000];
    szdate=size(datrange,2);
    mzFe_=zeros(szdate,1);mzGd_=zeros(szdate,1);
    angleFeGd_=zeros(szdate,1);
    Freqq_=zeros(szdate,1);
    for ctdat=1:szdate
        datname=sprintf('final%d.mat',datrange(ctdat));
        datname
        load(datname);
        if (0)
            figure;
            plot(t*1e12,mmAx,'linewidth',2)
            xlabel('time(ps)');ylabel('m_{Ax}')
        end
        if (1)
        subplot(2,1,1);
        %plot(tt*1e12,mmz(:,1),'linewidth',2);
        plot(t*1e12,mmAx,t*1e12,mmAy,t*1e12,mmAz,'linewidth',2);
        legend('mx','my','mz')
        xlabel('time(ps)');ylabel('m_A');
        subplot(2,1,2);
        %plot(tt*1e12,mmz(:,2),'linewidth',2);
        plot(t*1e12,mmBx,t*1e12,mmBy,t*1e12,mmBz,'linewidth',2);
        legend('mx','my','mz')
        xlabel('time(ps)');ylabel('m_B');
        end
        if (0)%FFT
            plotfft=0;
            nT0=size(mmAz,1);
            runTime0=runtime*1e9;%[ns]
            y=mmAz(:,1)';
            rminit=0.4;
            rmlast=0.01;
            Freqq_(ctdat)=FFT_module(nT0,runTime0,y,rminit,rmlast,plotfft);
        end
        close all;
        mzFe_(ctdat)=mmAz(end-1);
        mzGd_(ctdat)=mmBz(end-1);
        u=[mmAx(end-1),mmAy(end-1),mmAz(end-1)];
        v=[mmBx(end-1),mmBy(end-1),mmBz(end-1)];
        CosTheta = dot(u,v)/(norm(u)*norm(v));
        angleFeGd_(ctdat) = acosd(CosTheta);
    end
end

if (1)%two level
    datrang1=[-165:-5:-185];%Jc
    datrang2=[0:3:30];%alpNL
    sz1=size(datrang1,2);
    sz2=size(datrang2,2);
    angleFe_=zeros(sz1,sz2);
    for ctdat1=1:size(datrang1,2)
        for ctdat2=1:size(datrang2,2)
            datname=sprintf('final%d_%d',datrang1(ctdat1),datrang2(ctdat2));
            load(datname);
            [datrang1(ctdat1),datrang2(ctdat2)]
        if (0)
        subplot(2,1,1);
        %plot(tt*1e12,mmz(:,1),'linewidth',2);
        plot(t*1e12,mmAx,t*1e12,mmAy,t*1e12,mmAz,'linewidth',2);
        legend('mx','my','mz')
        xlabel('time(ps)');ylabel('m_A');
        subplot(2,1,2);
        %plot(tt*1e12,mmz(:,2),'linewidth',2);
        plot(t*1e12,mmBx,t*1e12,mmBy,t*1e12,mmBz,'linewidth',2);
        legend('mx','my','mz')
        xlabel('time(ps)');ylabel('m_B');
        end
        if (0)
        figure;
        plot(t*1e12,mmAx,'linewidth',2);
        xlabel('time(ps)');ylabel('m_A');
        end
        if(1)%terminal angle
            angleFe_(ctdat1,ctdat2)=acosd(mean(mmAx(end*0.9:end-1)));
        end
            close all
        end
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
if (0)
    Hex_Fe=6.5355;%[T]3meV->3meV/7.9293\mu_B=6.5355
    Hk_Fe_z=0.1742;
    Hk_Fe_x=23*0.1742;
    mub=9.274e-24;%[J/T]bohr magneton
    hbar=6.58211951440e-16;%eV.s
    ele=1.60200000000000e-19;%[C]
    gTM=2.2;%g-factor
    gamTM=gTM*mub/(hbar*ele);%1/(s.T)refer to "PRL 97, 217202 (2006), Jiang Xin"
    d=0.4e-9;%[m],lattice constant
    tz=d;
    musTM=2.217*mub;
    msTM=musTM/d^3;
    etaSTT=0.8;%spin hall angle
    JSTT_=[800:100:1000]*1e9;
    BDSTTTM=hbar/2*etaSTT*JSTT_/(msTM*tz);
    wSTT=gamTM*BDSTTTM*1e-12;%[THz]
    Freqnum=[0.0322033898305085;0.0440677966101695;0.0546610169491526]*14;%THz
    wperp=gamTM*Hk_Fe_x*1e-12;
    wpara=gamTM*Hk_Fe_z*1e-12;%[THz];
    alph=0.01;
    wE=gamTM*Hex_Fe*1e-12;%[THz]
    ws=linspace(0,1,100)*wperp;
    wn=wE*(1i*alph+sqrt((wperp+2*wpara-sqrt(wperp^2-4*ws.^2))./wE-alph^2));%acoustic mode
    figure;hold on
    plot(wSTT,Freqnum,'-*b')
    %plot(ws/wperp,real(wn)/2.4,'-r')
    plot(ws/wperp,real(wn),'-r')
    xlim([0 1]);%ylim([0 1.2])
    %legend()
    xlabel('{\omega}_s/{\omega}_{\perp}');ylabel('{\omega}(2{\pi}THz)')
end