function [ R,T ] = ScPnPv( vs,ps)
n=size(vs,2);
for i=1:n
    vs(:,i)=vs(:,i)/norm(vs(:,i));
end
%%
global Rs;
global Ts;
global Got;
global A;
global B;
global v;

Rs=zeros(3,3,200);
Ts=zeros(3,1,200);
Got=0;

as=zeros(1,10);
bs=zeros(1,10);
cs=zeros(1,10);
ds=zeros(1,10);

v=zeros(3,1);
V=zeros(3,3);

for i=1:n
    vi=vs(:,i);
    pi=ps(:,i);
    
    vx=vs(1,i);
    vy=vs(2,i);
    vz=vs(3,i);
    
    px=ps(1,i);
    py=ps(2,i);
    pz=ps(3,i);
    
    c_=[ px*vx^2 + py*vx*vy + pz*vx*vz, 2*py*vx*vz - 2*pz*vx*vy, 2*pz*vx^2 - 2*px*vz*vx, - 2*py*vx^2 + 2*px*vy*vx, px*vx^2 - py*vx*vy - pz*vx*vz, 2*py*vx^2 + 2*px*vy*vx, 2*pz*vx^2 + 2*px*vz*vx, py*vx*vy - px*vx^2 - pz*vx*vz, 2*py*vx*vz + 2*pz*vx*vy, pz*vx*vz - py*vx*vy - px*vx^2];
    as=as-c_;
    c_=[ py*vy^2 + px*vx*vy + pz*vy*vz, - 2*pz*vy^2 + 2*py*vz*vy, 2*pz*vx*vy - 2*px*vy*vz, 2*px*vy^2 - 2*py*vx*vy, px*vx*vy - py*vy^2 - pz*vy*vz, 2*px*vy^2 + 2*py*vx*vy, 2*px*vy*vz + 2*pz*vx*vy, py*vy^2 - px*vx*vy - pz*vy*vz, 2*pz*vy^2 + 2*py*vz*vy, pz*vy*vz - px*vx*vy - py*vy^2];
    bs=bs-c_;
    c_=[ pz*vz^2 + px*vx*vz + py*vy*vz, 2*py*vz^2 - 2*pz*vy*vz, 2*pz*vx*vz - 2*px*vz^2, 2*px*vy*vz - 2*py*vx*vz, px*vx*vz - pz*vz^2 - py*vy*vz, 2*px*vy*vz + 2*py*vx*vz, 2*px*vz^2 + 2*pz*vx*vz, py*vy*vz - px*vx*vz - pz*vz^2, 2*py*vz^2 + 2*pz*vy*vz, pz*vz^2 - px*vx*vz - py*vy*vz];
    cs=cs-c_;
    c_=[ px*vx + py*vy + pz*vz, 2*py*vz - 2*pz*vy, 2*pz*vx - 2*px*vz, 2*px*vy - 2*py*vx, px*vx - py*vy - pz*vz, 2*px*vy + 2*py*vx, 2*px*vz + 2*pz*vx, py*vy - px*vx - pz*vz, 2*py*vz + 2*pz*vy, pz*vz - py*vy - px*vx];
    ds=ds+c_;
    
    v=v+vs(:,i);
    V=V+eye(3)-vs(:,i)*vs(:,i)';    
end
C=eye(3)-v*v'/(v'*v);
A=(-V'*C*V-v*v')^(-1);
B=V'*C;

cost=getCost(vs,ps,as,bs,cs,ds,v,A,B);

e=-1;
for k=1:4
    %k
    [R_,T_,qz_,e_]=solveX1(cost(k,:),vs,ps,k);
    if e_<0
        continue;
    end
    if e<0 || e>e_
        R=R_;
        T=T_;
        e=e_;
    end
end
if e<0
    fprintf('ScPnP failed\r');
end

if n<=6
    R=Rs(:,:,1:Got);
    T=Ts(:,:,1:Got);
end
end

function [R,T,Q,e]=solveX1(cost,vs,ps,ind)
cost=cost/norm(cost);
R=eye(3);
T=zeros(3,1);
Q=zeros(4,1);

A1=xMqz1(cost);
A2=xMqz2(cost);
A3=xMqz3(cost);
A4=xMqz4(cost);
A5=xMqz5(cost);
A6=xMqz6(cost);
A7=xMqz7(cost);
A8=xMqz8(cost);
A9=xMqz9(cost);

Mn=17;
Ainv=A9^(-1);
B1=Ainv*A1;
B2=Ainv*A2;
B3=Ainv*A3;
B4=Ainv*A4;
B5=Ainv*A5;
B6=Ainv*A6;
B7=Ainv*A7;
B8=Ainv*A8;

I=eye(Mn);
O=zeros(Mn,Mn);
L=[ O,I,O,O,O,O,O,O;...
    O,O,I,O,O,O,O,O;...
    O,O,O,I,O,O,O,O;...
    O,O,O,O,I,O,O,O;...
    O,O,O,O,O,I,O,O;...
    O,O,O,O,O,O,I,O;...
    O,O,O,O,O,O,O,I;...
    -B1,-B2,-B3,-B4,-B5,-B6,-B7,-B8];

exclude=[1:15,17,18:25,27:30,32,34,35,36,38:40,42,44:46,49,52,55,56,59,62];

index=setdiff(1:136,exclude);
[v,d]=eig(L(index,index));
%[v,d]=eig(L);

e=-1;
for i=1:size(d,1)
    if abs(real(d(i,i)))<0.9
        continue;
    end
    if abs(imag(d(i,i)))>1E-8
        continue;
    end
    x1=1/real(d(i,i));
    Dx1 = A1*x1^8 +A2*x1^7 +A3*x1^6 +A4*x1^5 +A5*x1^4 +A6*x1^3 +A7*x1^2 +A8*x1^1 +A9;
    [u,s,v]=svd(Dx1); 
    x4= v(17,17);
    x3= v(15,17)/x4;
    x2= v(10,17)/x4;
    qw=1;
    qx=x1;
    qy=x2;
    qz=x3;
    if (ind==2)
       qw=x1;
       qx=1;
       qy=x2;
       qz=x3;
    elseif (ind ==3)
       qw=x1;
       qx=x2;
       qy=1;
       qz=x3;
    elseif(ind ==4)
       qw=x1;
       qx=x2;
       qy=x3;
       qz=1;
    end
    
    [R_,T_,qz_,e_]=recover(vs,ps,qw,qx,qy,qz);
    if e<0 || e>e_
        R=R_;
        T=T_;
        e=e_;
        Q=qz_;
    end
end
end
function [R,T,Q,e]=recover(vs,ps,qw,qx,qy,qz)
global Rs;
global Ts;
global Got;
global A;
global B;
global v;
n=size(vs,2);

q=[qw,qx,qy,qz];
Q=q/norm(q);
R=Q2R(q);

f1=zeros(3,1);
f2=0;
for i=1:n
    vi=vs(:,i);
    pi=ps(:,i);
    f1=f1-vi*vi'*R*pi;
    f2=f2+vi'*R*pi;
end

a=0;
b=0;
for i=1:n
    vi=vs(:,i);
    pi=ps(:,i);
    C_=eye(3)-vi*vi';
    A_=C_*(R*pi+A*B*f1+A*v*f2);
    a=a+A_'*A_;
    a_=C_*A*v;
    b=b+A_'*a_;
end
x4=a/b;
if x4<0
    x4=-x4;
end

T=A*(B*f1+v*f2-v*x4);

R=R/(q*q');
T=T/(q*q');

t=0;
a=0;
for i=1:n
    vi=vs(:,i);
    pi=ps(:,i);
    si=vi'*(R*pi+T);
    t=t+(vi-R*pi/si)/si;
    a=a+1/(si*si);
end

T=t/a;

Got=Got+1;
Rs(:,:,Got)=R;
Ts(:,:,Got)=T;
%
e=0;
for i=1:n
    vi=vs(:,i);
    vi=vi/norm(vi);
    pi=ps(:,i);
    si=vi'*(R*pi+T);
    temp=(R*pi+T)/si-vi;
    e=e+temp'*temp;
end
end



