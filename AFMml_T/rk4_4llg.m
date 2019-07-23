%% 4th order Runge Kutta method, for LLG calaulation in both IMA, PMA MTJ with FLT and DLT
% usage: add path which contain this file, call the function
% don't create the same function in new project

% call this function using: [,]=rk4_4llg(,)

%zzf,March.18,19.2016;
%1.changed based on PMA;2.add in FL torque
%% input
% Demag_, 3 by 3 matrix
% tstep is time step, unit [s]
% totstep is total number of steps
% m_init is initial magnetization, it is a 1-by-3 matrix, unit vector
% Ms: saturation magnetization, unit [emu/cm3]
% Hk: uniaxial anisotropy field, one value unit [tesla]
% Hext: applied field, 1-by-3 vector, unit [tesla]
% alp: damping constant
% P: polarization of FL and PL, currently only support same for both layer

% psj: unit 1-by-3 vector, spin flux polarization,
% note in STT the reflection type is opposite to m_pin_layer

% dimension FL_length,FL_width,FL_thickness, unit [nm]

%% output
%mmx,mmy,mmz: magnetization component, unit vector
%tt: simulation time list, unit [ns]
%Icri: critical current for switching unit:[Ampere]
if dimensionlessLLG
    Hk_=Hk;
    Hk=[1*(FL_width<FL_length)*Hk,1*(FL_width>FL_length)*Hk,0];
    tau_c=(g*Hk(2*(FL_width>FL_length)+1*(FL_width<FL_length)))/(1+alp^2); %natural time constant 1/s
    scal=1;
else
    %Hk=[1,1,1];%normalization purpose
    tau_c=1;%[1/s]
    scal=1;%scale parameter
end
ts1=tstep*tau_c; %time step

ct1=1; %count 1
t=linspace(0,runtime,totstep);
mmx=zeros(totstep,1);mmy=zeros(totstep,1);mmz=zeros(totstep,1);
Lx=zeros(totstep,1);Ly=zeros(totstep,1);Lz=zeros(totstep,1);
mmAx=zeros(totstep,1);mmAy=zeros(totstep,1);mmAz=zeros(totstep,1);
mmBx=zeros(totstep,1);mmBy=zeros(totstep,1);mmBz=zeros(totstep,1);

mmx(1,1)=m_init1(1)+m_init2(1);%m=(mA+mB)/2 [1]
mmy(1,1)=m_init1(2)+m_init2(2);
mmz(1,1)=m_init1(3)+m_init2(3);
Lx(1,1)=m_init1(1)-m_init2(1);%m=(mA-mB)/2 [1]
Ly(1,1)=m_init1(2)-m_init2(2);
Lz(1,1)=m_init1(3)-m_init2(3);
mmAx(1,1)=m_init1(1);%first sublattice
mmAy(1,1)=m_init1(2);
mmAz(1,1)=m_init1(3);
mmBx(1,1)=m_init2(1);%second sublattice
mmBy(1,1)=m_init2(2);
mmBz(1,1)=m_init2(3);
while ct1<totstep
    mm1=[mmx(ct1,1),mmy(ct1,1),mmz(ct1,1)];
    LL1=[Lx(ct1,1),Ly(ct1,1),Lz(ct1,1)];
    mmA1=[mmAx(ct1,1),mmAy(ct1,1),mmAz(ct1,1)];
    mmB1=[mmBx(ct1,1),mmBy(ct1,1),mmBz(ct1,1)];
    
    mmm=mm1;
    L=LL1;
    mmA=(mmm+L)/2;
    mmB=(mmm-L)/2;
    [~,sttdlt,~,~,~]=field_eta(mmm,Hk,Demag_,Hext,jc_STT,...
        tFL,Ms,facFLT_SHE,K12Dipole,mmmPL,PolFL,LFL,WFL,facFLT_STT,...
        thetaSH,tHM,lambdaSF,jc_SOT,TT,alp,tstep,thermalnois);
    fL=(2*wE*L+[w_para*L(1),0,0]-[0,0,w_perp*L(3)]);
    fm=([w_para*mmm(1),0,0]-[0,0,w_perp*mmm(3)]);
    ws=sttdlt*PolSTT;
    if ct1==1
        mmm_Fe=m_init1;
        mmm_Gd=m_init2;
        Hk_Fe_z=-w_perp;
        Hk_Gd_z=Hk_Fe_z;
        Hk_Fe_x=w_para;
        Hk_Gd_x=Hk_Fe_x;
        Hex_Fe=wE;
        Hex_Gd=Hex_Fe;
        Ms_Fe=Ms;
        Ms_Gd=Ms_Fe;
        [hh_Fe,sttdlt_Fe,sttflt_Fe,sotdlt_Fe,sotflt_Fe,hh_Gd,sttdlt_Gd,sttflt_Gd...
            ,sotdlt_Gd,sotflt_Gd]=field_eta_m1m2(mmm_Fe,mmm_Gd,Hk_Fe_z,Hk_Gd_z,Hk_Fe_x,Hk_Gd_x,Demag_,Hext,jc_STT,...
            tFL,Ms_Fe,Ms_Gd,facFLT_SHE,K12Dipole,mmmPL,PolFL,LFL,WFL,facFLT_STT,...
            thetaSH,tHM,lambdaSF,jc_SOT,TT,alp,tstep,thermalnois,Hex_Fe,Hex_Gd);
        dmdt_Fe=LLG_solver_m1m2(alp,mmm_Fe,hh_Fe,polSOT,PolSTT,sttdlt_Fe,sttflt_Fe,sotdlt_Fe,sotflt_Fe);
        dmdt_Gd=LLG_solver_m1m2(alp,mmm_Gd,hh_Gd,polSOT,PolSTT,sttdlt_Gd,sttflt_Gd,sotdlt_Gd,sotflt_Gd);
        scaltmp=gam/(1+alp^2);
        dLdt_tmp=scaltmp*(dmdt_Fe-dmdt_Gd);
        dmdt_tmp=scaltmp*(dmdt_Fe+dmdt_Gd);
        dLdt=dLdt_tmp;
        dmdt=dmdt_tmp;
        kk1m=scal*dmdt;
        kk1L=scal*dLdt;
    else
        dLdt_tmp=dLdt;
        dmdt_tmp=dmdt;
        [dmdt,dLdt]=LLG_solver(alp,mmm,L,fL,fm,dLdt_tmp,dmdt_tmp,ws,gam);
        kk1m=scal*dmdt;
        kk1L=scal*dLdt;
    end

    mmm=mm1+kk1m*ts1/2;
    L=LL1+kk1L*ts1/2;
    mmA=(mmm+L)/2;
    mmB=(mmm-L)/2;
    [~,sttdlt,~,~,~]=field_eta(mmm,Hk,Demag_,Hext,jc_STT,...
        tFL,Ms,facFLT_SHE,K12Dipole,mmmPL,PolFL,LFL,WFL,facFLT_STT,...
        thetaSH,tHM,lambdaSF,jc_SOT,TT,alp,tstep,thermalnois);
    fL=(2*wE*L+[w_para*L(1),0,0]-[0,0,w_perp*L(3)]);
    fm=([w_para*mmm(1),0,0]-[0,0,w_perp*mmm(3)]);
    ws=sttdlt*PolSTT;
    dLdt_tmp=dLdt;
    dmdt_tmp=dmdt;
    [dmdt,dLdt]=LLG_solver(alp,mmm,L,fL,fm,dLdt_tmp,dmdt_tmp,ws,gam);
    kk2m=scal*dmdt;
    kk2L=scal*dLdt;
    
    mmm=mm1+kk2m*ts1/2;
    L=LL1+kk2L*ts1/2;
    mmA=(mmm+L)/2;
    mmB=(mmm-L)/2;
    [~,sttdlt,~,~,~]=field_eta(mmm,Hk,Demag_,Hext,jc_STT,...
        tFL,Ms,facFLT_SHE,K12Dipole,mmmPL,PolFL,LFL,WFL,facFLT_STT,...
        thetaSH,tHM,lambdaSF,jc_SOT,TT,alp,tstep,thermalnois);
    fL=(2*wE*L+[w_para*L(1),0,0]-[0,0,w_perp*L(3)]);
    fm=([w_para*mmm(1),0,0]-[0,0,w_perp*mmm(3)]);
    ws=sttdlt*PolSTT;
    dLdt_tmp=dLdt;
    dmdt_tmp=dmdt;
    [dmdt,dLdt]=LLG_solver(alp,mmm,L,fL,fm,dLdt_tmp,dmdt_tmp,ws,gam);
    kk3m=scal*dmdt;
    kk3L=scal*dLdt;
    
    mmm=mm1+kk3m*ts1;
    L=LL1+kk3L*ts1;
    mmA=(mmm+L)/2;
    mmB=(mmm-L)/2;
    [~,sttdlt,~,~,~]=field_eta(mmm,Hk,Demag_,Hext,jc_STT,...
        tFL,Ms,facFLT_SHE,K12Dipole,mmmPL,PolFL,LFL,WFL,facFLT_STT,...
        thetaSH,tHM,lambdaSF,jc_SOT,TT,alp,tstep,thermalnois);
    fL=(2*wE*L+[w_para*L(1),0,0]-[0,0,w_perp*L(3)]);
    fm=([w_para*mmm(1),0,0]-[0,0,w_perp*mmm(3)]);
    ws=sttdlt*PolSTT;
    dLdt_tmp=dLdt;
    dmdt_tmp=dmdt;
    [dmdt,dLdt]=LLG_solver(alp,mmm,L,fL,fm,dLdt_tmp,dmdt_tmp,ws,gam);
    kk4m=scal*dmdt;
    kk4L=scal*dLdt;
    
    mn1=mm1+ts1/6*(kk1m+2*kk2m+2*kk3m+kk4m);
    Ln1=LL1+ts1/6*(kk1L+2*kk2L+2*kk3L+kk4L);
    
    mmAx_tmp=(mn1(1)+Ln1(1))/2;
    mmAy_tmp=(mn1(2)+Ln1(2))/2;
    mmAz_tmp=(mn1(3)+Ln1(3))/2;
    mmBx_tmp=(mn1(1)-Ln1(1))/2;
    mmBy_tmp=(mn1(2)-Ln1(2))/2;
    mmBz_tmp=(mn1(3)-Ln1(3))/2;
    mmA_tmp2=[mmAx_tmp,mmAy_tmp,mmAz_tmp];
    mmB_tmp2=[mmBx_tmp,mmBy_tmp,mmBz_tmp];
    mmA_tmp3=mmA_tmp2/norm(mmA_tmp2);
    mmB_tmp3=mmB_tmp2/norm(mmB_tmp2);
    mmAx(ct1+1,1)=mmA_tmp3(1);
    mmAy(ct1+1,1)=mmA_tmp3(2);
    mmAz(ct1+1,1)=mmA_tmp3(3);
    mmBx(ct1+1,1)=mmB_tmp3(1);
    mmBy(ct1+1,1)=mmB_tmp3(2);
    mmBz(ct1+1,1)=mmB_tmp3(3);
    
    mmx(ct1+1,1)=mmAx(ct1+1,1)+mmBx(ct1+1,1);
    mmy(ct1+1,1)=mmAy(ct1+1,1)+mmBy(ct1+1,1);
    mmz(ct1+1,1)=mmAz(ct1+1,1)+mmBz(ct1+1,1);
    Lx(ct1+1,1)=mmAx(ct1+1,1)-mmBx(ct1+1,1);
    Ly(ct1+1,1)=mmAy(ct1+1,1)-mmBy(ct1+1,1);
    Lz(ct1+1,1)=mmAz(ct1+1,1)-mmBz(ct1+1,1);
    
    ct1=ct1+1;
end

Lx_tmp=Lx(1:savestep:end);
Ly_tmp=Ly(1:savestep:end);
Lz_tmp=Lz(1:savestep:end);
mmx_tmp=mmx(1:savestep:end);
mmy_tmp=mmy(1:savestep:end);
mmz_tmp=mmz(1:savestep:end);
mmAx_tmp=mmAx(1:savestep:end);
mmAy_tmp=mmAy(1:savestep:end);
mmAz_tmp=mmAz(1:savestep:end);
mmBx_tmp=mmBx(1:savestep:end);
mmBy_tmp=mmBy(1:savestep:end);
mmBz_tmp=mmBz(1:savestep:end);
t_tmp=t(1:savestep:end);
clear Lx Ly Lz mmx mmy mmz mmAx mmAy mmAz mmBx mmBy mmBz t
Lx=Lx_tmp;
Ly=Ly_tmp;
Lz=Lz_tmp;
mmx=mmx_tmp; 
mmy=mmy_tmp; 
mmz=mmz_tmp; 
mmAx=mmAx_tmp; 
mmAy=mmAy_tmp; 
mmAz=mmAz_tmp; 
mmBx=mmBx_tmp; 
mmBy=mmBy_tmp;  
mmBz=mmBz_tmp;  
t=t_tmp;
clear Lx_tmp Ly_tmp Lz_tmp mmx_tmp mmy_tmp mmz_tmp 
clear mmAx_tmp mmAy_tmp mmAz_tmp mmBx_tmp mmBy_tmp mmBz_tmp t_tmp


