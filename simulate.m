function [ps,vs,Ti,Ri] = simulate(s,n,az)
Ti=rand(3,1)*2-1;
qi=rand(4,1)*2-1;
qi=qi/norm(qi);
Ri=Q2R(qi);

ps=rand(3,n)*2-1;
pc=mean(ps,2);
for i=1:n
   ps(:,i)=ps(:,i)-pc;
end

%coplanar
% [uu,ss,vv]=svd(ps);
% ss(3,3)=0;
% ps=uu*ss*vv';
% pc=mean(ps,2);
% for i=1:n
%    ps(:,i)=ps(:,i)-pc;
% end
%coplanar


phz=inf;
for i=1:n
   ph=Ri*ps(:,i)+Ti;
   if phz>ph(3)
       phz=ph(3);
   end
end
if az>phz
    Ti(3)=Ti(3)+az-phz;% ScPnPv does not need this
end

vs=zeros(3,n);
for i=1:n
    ph=Ri*ps(:,i)+Ti;
    
    % for perspective camera
%     e=normrnd(0,s,2,1);
%     vs(:,i) = ph/ph(3)+[e;0];
    
    % for sphere camera
    e=normrnd(0,s,3,1);
    vs(:,i) = ph/norm(ph)+e;
    vs(:,i)=vs(:,i)/vs(3,i);
end
end