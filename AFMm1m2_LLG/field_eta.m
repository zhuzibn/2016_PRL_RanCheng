%% function for effective field calculation
% usage: add path which contain this file, call the function
% don't create the same function in new project 
%input
%1. mmm, magnetization, [1x3] vector, [emu/cm3]
%2. Hk, crystalline anisotropy field, double, [Tesla]
%3. Demag_, demagnetizing tensor, [3x3] matrix
%4. Hext, external field, [1x3] vector, [Tesla]
%5. jc_STT, spin current density, double, [A/m2]
%6. tFL, free layer thickness, double, [m]
%7. Ms, saturation magnetization, double, [emu/cm3]
%8. facFLT_SHE,ratio of FLT/DLT
%9. K12Dipole, dipole tensor, [3x3] matrix
%10. mmmPL: magnetization of PL
%11. PolFL:polarization of FL
%12. lFL: length of FL [m]
%13. wFL: width of FL [m]
%14. facFLT_STT:ratio of FLT/DLT in STT 
%15. thetaSH:spin hall angle
%16. tHM:[m] thickness of HM
%17. lambdaSF: [m]spin diffusion length
%18. jc_SOT:[A/m2] SOT current density
%19. TT:[K] Temperature
%20. alp:damping constant, dimensionless
%21. tstep: [s] time step
%22. thermalnois: flag for thermal noise
%output
%1. hh,total effective field (include SOT FLT), [1x3] vector, [Tesla]
%2. sttdlt, STT DLT cofficient, double, [Tesla]
%3. sttflt, STT DLT cofficient, double, [Tesla]
%4. sotdlt, STT DLT cofficient, double, [Tesla]
%5. sotflt, STT DLT cofficient, double, [Tesla]
function [hh_Fe,sttdlt_Fe,sttflt_Fe,sotdlt_Fe,sotflt_Fe,hh_Gd,sttdlt_Gd,sttflt_Gd,sotdlt_Gd,sotflt_Gd]=field_eta(mmm_Fe,mmm_Gd,Hk_Fe_z,Hk_Gd_z,Hk_Fe_x,Hk_Gd_x,Demag_,Hext,jc_STT,...
    tFL,Ms_Fe,Ms_Gd,facFLT_SHE,K12Dipole,mmmPL,PolFL,lFL,wFL,facFLT_STT,...
    thetaSH,tHM,lambdaSF,jc_SOT,TT,alp,tstep,thermalnois,Hex_Fe,Hex_Gd)
conf_file();%load configuration
constantfile();%load constant
hk_Fe=[Hk_Fe_x,0,Hk_Fe_z].*mmm_Fe;
hk_Gd=[Hk_Gd_x,0,Hk_Gd_z].*mmm_Gd;

Hd_Fe=4*pi*Ms_Fe*Demag_*1e-4;%Tesla
Hd_Gd=4*pi*Ms_Gd*Demag_*1e-4;%Tesla

hd_Fe=(-Hd_Fe*mmm_Fe')';%([3x3]*[3x1])'=[1x3]
hd_Gd=(-Hd_Gd*mmm_Gd')';%([3x3]*[3x1])'=[1x3]

hex_Fe=-Hex_Fe*mmm_Gd;
hex_Gd=-Hex_Gd*mmm_Fe;

hext=Hext;

if dipolee
    hdipole = 4*pi*1e-7*(K12Dipole*Ms_Fe*1e3*mmmPL')';%[Tesla]
else
    hdipole=[0,0,0];
end

Jp_Fe=2*tFL*(Ms_Fe*1e3)/hbar;
Jp_Gd=2*tFL*(Ms_Gd*1e3)/hbar;
if STT_DLT
    efficiencyselect=3;
    switch efficiencyselect%only for IMA, to modify to fit for PMA
        case 1
            b=1/(-4+((1+P)^3)*(3+mmm_Fe(2)*(wPL>lPL)+mmm_Fe(1)*(wPL<lPL)))/(4*(P^(3/2)));
            %Slonswski torque efficiency for GMR
        case 2%Slonswski torque efficiency for TMR
            b=PolFL/(1+(PolFL^2)*(dot(mmm_Fe,mmmPL)));
        case 3
            b=0.8; % fixed torque efficiency
    end
    sttdlt_Fe=jc_STT/Jp_Fe*b;
    sttdlt_Gd=jc_STT/Jp_Gd*b;
    sttflt_Fe=facFLT_STT*sttdlt_Fe;
    sttflt_Gd=facFLT_STT*sttdlt_Gd;
else
    sttdlt_Fe=0;
    sttdlt_Gd=0;
    sttflt_Fe=0;
    sttflt_Gd=0;
end
if SOT_DLT
    sotdlt_Fe=thetaSH*jc_SOT/Jp_Fe*(1-sech(tHM/lambdaSF));%to modify to auto get easy (y) axis
    sotdlt_Gd=thetaSH*jc_SOT/Jp_Gd*(1-sech(tHM/lambdaSF));
    sotflt_Fe=facFLT_SHE*sotdlt_Fe;
    sotflt_Gd=facFLT_SHE*sotdlt_Gd;
else
    sotdlt_Fe=0;
    sotdlt_Gd=0;
    sotflt_Fe=0;
    sotflt_Gd=0;
end
%% thermal fluctuation
if thermalnois==1%1(0) (not) enable thermal noise
    hthermtmp=sqrt(2*kb*TT*alp/(lFL*wFL*tFL*Ms_Fe*1e3*gam*(1+alp^2)*tstep));%[T]
    hthermx=normrnd(0,hthermtmp);hthermy=normrnd(0,hthermtmp);hthermz=normrnd(0,hthermtmp);
    htherm=[hthermx,hthermy,hthermz];
else
    htherm=[0,0,0];
end
hh_Fe=hk_Fe+hd_Fe+hext+hdipole+htherm+hex_Fe; %total field
hh_Gd=hk_Gd+hd_Gd+hext+hdipole+htherm+hex_Gd;
end