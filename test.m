clc;
clear;
format long
addpath ScPnPv;
addpath ScPnPz;

sigma=[0.01:0.01:0.08];
NN=size(sigma,2);
EN = 500;

Ed=zeros(2,NN);
Esc=zeros(2,NN);
for nn=1:NN
    esc=zeros(2,EN);
    ed=zeros(2,EN);
    for en=1:EN
        [ps,vs,Ti,Ri] = simulate(sigma(nn),12,1);
        %----------------------------------------------------------------------
        [ Rsc,Tsc ] = ScPnPv(vs,ps);
        esc(:,en) = ErrorM( Ri,Rsc,Ti,Tsc );
        [ Rd,Td ] = ScPnPz(vs,ps);
        ed(:,en) = ErrorM( Ri,Rd,Ti,Td );
        %----------------------------------------------------------------------
    end
    for k=1:2
        Esc(k,nn) = mean(esc(k,:));
        Ed(k,nn) = mean(ed(k,:));
    end
end