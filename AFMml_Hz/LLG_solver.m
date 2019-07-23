%% LLG equation with precession term, damping term, spin current
%refer to Eq. (4,5) in Gomonay's 2010 PRB paper for the complete equation
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
function [dmdt,dLdt]=LLG_solver(alp,mmm,L,fL,fm,dLdt_tmp,dmdt_tmp,ws,alp_NL)
% call this function by feval(@(t,m) LLG_solver(t,m,Hk,alpha),t0,m0)
% t0 is the initial value of t
% m0 is the initial value of m
dmdt=0.5*(cross(fL,L)+cross(fm,mmm))+alp/2*(cross(mmm,dmdt_tmp)+cross(L,dLdt_tmp))+...
    0.5*(cross(L,cross(ws,L))+cross(mmm,cross(ws,mmm))+...
    alp_NL*(L(3)^2*cross(L,dLdt_tmp)-dLdt_tmp(3)*cross([0,0,1],L)));
dLdt=0.5*(cross(fm,L)+cross(fL,mmm))+alp/2*(cross(mmm,dLdt_tmp)+cross(L,dmdt_tmp))+...
    0.5*(cross(L,cross(ws,mmm))+cross(mmm,cross(ws,L)));
end