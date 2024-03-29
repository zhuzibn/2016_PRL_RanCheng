%% LLG equation with precession term, damping term, spin current
% usage: add path which contain this file, call the function
% don't create the same function in new project 
%1.alp:damping constant,value
%2.mmm:magnetization, 1-by-3 matrix
%3.hh:effective field, 1-by-3 matrix
%4.pSOT:SOT polarization, 1-by-3 matrix
%5.pSTT:STT polarization, 1-by-3 matrix
%6.sttdlt:strength of STT DLT,value
%7.sttflt:strength of STT DLT,value
%8.sotdlt:strength of STT DLT,value
%9.sotflt:strength of STT DLT,value
function dmdt=LLG_solver(gam,alp,mmm,hh,pSTT,sttdlt,sttflt,dmdt_tmp)
% call this function by feval(@(t,m) LLG_solver(t,m,Hk,alpha),t0,m0)
% t0 is the initial value of t
% m0 is the initial value of m

dmdt=-gam*cross(mmm,hh)+alp*cross(mmm,dmdt_tmp)-...
    gam*sttdlt*cross(mmm,cross(mmm,pSTT))+gam*sttflt*cross(mmm,pSTT);
    
end