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
    tau_c=1;
    scal=1;%scale parameter
end
ts1=tstep*tau_c; %time step

ct1=1; %count 1
t=linspace(0,runtime,totstep);
mmx=zeros(totstep,2);%(:,1)is Fe, (:,2)is Gd
mmy=zeros(totstep,2);
mmz=zeros(totstep,2);
mmx(1,1)=ma_init(1);mmy(1,1)=ma_init(2);mmz(1,1)=ma_init(3);
mmx(1,2)=mb_init(1);mmy(1,2)=mb_init(2);mmz(1,2)=mb_init(3);
while ct1<totstep
    mmFe=[mmx(ct1,1),mmy(ct1,1),mmz(ct1,1)]; %Fe
    mmGd=[mmx(ct1,2),mmy(ct1,2),mmz(ct1,2)]; %Gd
    %% current calc
    if (0) %unit conversion Tesla-->A/m,debug use
        hhtmp=(hh*Hk(2)*1e7)/4/pi;
        hk=(hk*Hk(2)*1e7)/4/pi;
        hd=(hd*Hk(2)*1e7)/4/pi;
        hext=(hext*Hk(2)*1e7)/4/pi;
        hdipole=(hdipole*Hk(2)*1e7)/4/pi;
    end
    %% unit convension:
    %e_tmp:[e],unit electron charge
    %hbar_tmp:[ev.s]
    %d_tmp:[m],FL thickness
    %Hk_tmp:[Gauss]
    mmm_Fe=mmFe;
    mmm_Gd=mmGd;
    [hh_Fe,sttdlt_Fe,sttflt_Fe,sotdlt_Fe,sotflt_Fe,hh_Gd,sttdlt_Gd,sttflt_Gd...
        ,sotdlt_Gd,sotflt_Gd]=field_eta(mmm_Fe,mmm_Gd,Hk_Fe_z,Hk_Gd_z,Hk_Fe_x,Hk_Gd_x,Demag_,Hext,jc_STT,...
        tFL,Ms_Fe,Ms_Gd,facFLT_SHE,K12Dipole,mmmPL,PolFL,LFL,WFL,facFLT_STT,...
        thetaSH,tHM,lambdaSF,jc_SOT,TT,alp,tstep,thermalnois,Hex_Fe,Hex_Gd);
    if ct1==1
    dmdt_Fe=gam/(1+alp^2)*LLG_solver_tmp(alp,mmm_Fe,hh_Fe,polSOT,PolSTT,sttdlt_Fe,sttflt_Fe,sotdlt_Fe,sotflt_Fe);
    dmdt_Gd=gam/(1+alp^2)*LLG_solver_tmp(alp,mmm_Gd,hh_Gd,polSOT,PolSTT,sttdlt_Gd,sttflt_Gd,sotdlt_Gd,sotflt_Gd);        
    else
    dmdt_Fe_tmp=dmdt_Fe;    
    dmdt_Gd_tmp=dmdt_Gd;
    dmdt_Fe=LLG_solver(gam,alp,mmm_Fe,hh_Fe,PolSTT,sttdlt_Fe,sttflt_Fe,dmdt_Fe_tmp);  
    dmdt_Gd=LLG_solver(gam,alp,mmm_Gd,hh_Gd,PolSTT,sttdlt_Gd,sttflt_Gd,dmdt_Gd_tmp);     
    end
    kk1_Fe=scal*dmdt_Fe;
    kk1_Gd=scal*dmdt_Gd;
    
    mmm_Fe=mmFe+kk1_Fe*ts1/2;
    mmm_Gd=mmGd+kk1_Gd*ts1/2;
    [hh_Fe,sttdlt_Fe,sttflt_Fe,~,~,hh_Gd,sttdlt_Gd,sttflt_Gd...
        ,~,~]=field_eta(mmm_Fe,mmm_Gd,Hk_Fe_z,Hk_Gd_z,Hk_Fe_x,Hk_Gd_x,Demag_,Hext,jc_STT,...
        tFL,Ms_Fe,Ms_Gd,facFLT_SHE,K12Dipole,mmmPL,PolFL,LFL,WFL,facFLT_STT,...
        thetaSH,tHM,lambdaSF,jc_SOT,TT,alp,tstep,thermalnois,Hex_Fe,Hex_Gd);
    dmdt_Fe_tmp=dmdt_Fe;    
    dmdt_Gd_tmp=dmdt_Gd;
    dmdt_Fe=LLG_solver(gam,alp,mmm_Fe,hh_Fe,PolSTT,sttdlt_Fe,sttflt_Fe,dmdt_Fe_tmp);  
    dmdt_Gd=LLG_solver(gam,alp,mmm_Gd,hh_Gd,PolSTT,sttdlt_Gd,sttflt_Gd,dmdt_Gd_tmp);     
    kk2_Fe=scal*dmdt_Fe;
    kk2_Gd=scal*dmdt_Gd;
    
    mmm_Fe=mmFe+kk2_Fe*ts1/2;
    mmm_Gd=mmGd+kk2_Gd*ts1/2;
    [hh_Fe,sttdlt_Fe,sttflt_Fe,~,~,hh_Gd,sttdlt_Gd,sttflt_Gd...
        ,~,~]=field_eta(mmm_Fe,mmm_Gd,Hk_Fe_z,Hk_Gd_z,Hk_Fe_x,Hk_Gd_x,Demag_,Hext,jc_STT,...
        tFL,Ms_Fe,Ms_Gd,facFLT_SHE,K12Dipole,mmmPL,PolFL,LFL,WFL,facFLT_STT,...
        thetaSH,tHM,lambdaSF,jc_SOT,TT,alp,tstep,thermalnois,Hex_Fe,Hex_Gd);
    dmdt_Fe_tmp=dmdt_Fe;    
    dmdt_Gd_tmp=dmdt_Gd;
    dmdt_Fe=LLG_solver(gam,alp,mmm_Fe,hh_Fe,PolSTT,sttdlt_Fe,sttflt_Fe,dmdt_Fe_tmp);  
    dmdt_Gd=LLG_solver(gam,alp,mmm_Gd,hh_Gd,PolSTT,sttdlt_Gd,sttflt_Gd,dmdt_Gd_tmp);     
    kk3_Fe=scal*dmdt_Fe;
    kk3_Gd=scal*dmdt_Gd;
    
    mmm_Fe=mmFe+kk3_Fe*ts1;
    mmm_Gd=mmGd+kk3_Gd*ts1;
    [hh_Fe,sttdlt_Fe,sttflt_Fe,~,~,hh_Gd,sttdlt_Gd,sttflt_Gd...
        ,~,~]=field_eta(mmm_Fe,mmm_Gd,Hk_Fe_z,Hk_Gd_z,Hk_Fe_x,Hk_Gd_x,Demag_,Hext,jc_STT,...
        tFL,Ms_Fe,Ms_Gd,facFLT_SHE,K12Dipole,mmmPL,PolFL,LFL,WFL,facFLT_STT,...
        thetaSH,tHM,lambdaSF,jc_SOT,TT,alp,tstep,thermalnois,Hex_Fe,Hex_Gd);
    dmdt_Fe_tmp=dmdt_Fe;    
    dmdt_Gd_tmp=dmdt_Gd;
    dmdt_Fe=LLG_solver(gam,alp,mmm_Fe,hh_Fe,PolSTT,sttdlt_Fe,sttflt_Fe,dmdt_Fe_tmp);  
    dmdt_Gd=LLG_solver(gam,alp,mmm_Gd,hh_Gd,PolSTT,sttdlt_Gd,sttflt_Gd,dmdt_Gd_tmp);     
    kk4_Fe=scal*dmdt_Fe;
    kk4_Gd=scal*dmdt_Gd;
    
    mnFe=mmFe+ts1/6*(kk1_Fe+2*kk2_Fe+2*kk3_Fe+kk4_Fe);
    mnGd=mmGd+ts1/6*(kk1_Gd+2*kk2_Gd+2*kk3_Gd+kk4_Gd);
    mnFe=mnFe/norm(mnFe);
    mnGd=mnGd/norm(mnGd);
    mmx(ct1+1,1)=mnFe(1);mmy(ct1+1,1)=mnFe(2);mmz(ct1+1,1)=mnFe(3);
    mmx(ct1+1,2)=mnGd(1);mmy(ct1+1,2)=mnGd(2);mmz(ct1+1,2)=mnGd(3);
    
    ct1=ct1+1;
end
if dimensionlessLLG
    tt=t/tau_c*1e9;%unit[ns]
else
    tt=t;
end
mmx_tmp=mmx(1:savestep:end,:);
mmy_tmp=mmy(1:savestep:end,:);
mmz_tmp=mmz(1:savestep:end,:);
tt_tmp=tt(1:savestep:end);
clear mmx mmy mmz tt t
mmx=mmx_tmp;mmy=mmy_tmp;mmz=mmz_tmp;tt=tt_tmp;
clear mmx_tmp mmy_tmp mmz_tmp tt_tmp