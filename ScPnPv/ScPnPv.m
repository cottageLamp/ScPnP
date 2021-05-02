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
    [R_,T_,qz_,e_]=solveQz(cost(k,:),vs,ps,setdiff(1:4,k));
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

function [R,T,Q,e]=solveQz(cost,vs,ps,ind)
cost=cost/norm(cost);
R=eye(3);
T=zeros(3,1);
Q=zeros(4,1);

A1=zMqz1(cost);
A2=zMqz2(cost);
A3=zMqz3(cost);
A4=zMqz4(cost);
A5=zMqz5(cost);
A6=zMqz6(cost);
A7=zMqz7(cost);
A8=zMqz8(cost);
A9=zMqz9(cost);
A10=zMqz10(cost);

Mn=22;
Ainv=A10^(-1);
B1=Ainv*A1;
B2=Ainv*A2;
B3=Ainv*A3;
B4=Ainv*A4;
B5=Ainv*A5;
B6=Ainv*A6;
B7=Ainv*A7;
B8=Ainv*A8;
B9=Ainv*A9;
I=eye(Mn);
O=zeros(Mn,Mn);
L=[O,I,O,O,O,O,O,O,O;...
    O,O,I,O,O,O,O,O,O;...
    O,O,O,I,O,O,O,O,O;...
    O,O,O,O,I,O,O,O,O;...
    O,O,O,O,O,I,O,O,O;...
    O,O,O,O,O,O,I,O,O;...
    O,O,O,O,O,O,O,I,O;...
    O,O,O,O,O,O,O,O,I;...
    -B1,-B2,-B3,-B4,-B5,-B6,-B7,-B8,-B9];

exclude=[1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	21	22	23	24	25	26	27	28	29	30	32	33	34	35	36	37	38	40	41	43	44	45	46	48	49	50	52	54	55	56	57	58	60	62	63	66	70	71	74	77	78	80	82	85	88	102	107];
%temp=L(177:198,exclude);
%if norm(temp)>0.0000001
 %   printf('error\r');
%end
%[~,d]=eig(L);
index=setdiff(1:198,exclude);
[v,d]=eig(L(index,index));

e=-1;
for i=1:size(d,1)
    if abs(real(d(i,i)))<0.9
        continue;
    end
    if abs(imag(d(i,i)))>1E-8
        continue;
    end
    qz=1/real(d(i,i));
    Dqz = A1*qz^9 +A2*qz^8 +A3*qz^7 +A4*qz^6 +A5*qz^5 +A6*qz^4 +A7*qz^3 +A8*qz^2 +A9*qz + A10;
    [u,s,v]=svd(Dqz); 
    x4= v(21,22);
    qy= v(18,22)/x4;
    qx= v(10,22)/x4;
    
    [R_,T_,qz_,e_]=recover(vs,ps,qz,qx,qy,ind);
    if e<0 || e>e_
        R=R_;
        T=T_;
        e=e_;
        Q=qz_;
    end
end
end
function [R,T,Q,e]=solveQx(cost,vs,ps,qz,ind)
R=eye(3);
T=zeros(3,1);
Q=zeros(4,1);
As=zMqx(cost,qz);
Mn=5;
A1=reshape(As(1,:,:),Mn,Mn);
A2=reshape(As(2,:,:),Mn,Mn);
A3=reshape(As(3,:,:),Mn,Mn);
A4=reshape(As(4,:,:),Mn,Mn);
A5=reshape(As(5,:,:),Mn,Mn);
A6=reshape(As(6,:,:),Mn,Mn);
A7=reshape(As(7,:,:),Mn,Mn);
A8=reshape(As(8,:,:),Mn,Mn);
%   A=A1*x^9+A2*x^8+A3*x^7+A4*x^6+A5*x^5+A6*x^4+A7*x^3+A8*x^2+A9*x^1+A10;
%   [uu,ss,vv]=svd(A);
%   ss(22,22)
%[uu,ss,vv]=svd(A8);
%[qz,ss(5,5)]

A8=A8^(-1);
A1=A8*A1;
A2=A8*A2;
A3=A8*A3;
A4=A8*A4;
A5=A8*A5;
A6=A8*A6;
A7=A8*A7;
I=eye(Mn);
O=zeros(Mn,Mn);
L=[O,I,O,O,O,O,O;...
    O,O,I,O,O,O,O;...
    O,O,O,I,O,O,O;...
    O,O,O,O,I,O,O;...
    O,O,O,O,O,I,O;...
    O,O,O,O,O,O,I;...
    -A1,-A2,-A3,-A4,-A5,-A6,-A7];
[~,d]=eig(L);
e=-1;
for i=1:size(d,1)
    if abs(real(d(i,i)))<0.9
        continue;
    end
    if abs(imag(d(i,i)))>1E-8
        continue;
    end
    qx=1/real(d(i,i));
    [R_,T_,qz_,e_]=solveQy(cost,vs,ps,qz,qx,ind);
    if e_<0
        continue;
    end
    if e<0 || e>e_
        R=R_;
        T=T_;
        e=e_;
        Q=qz_;
    end
end
end
function [R,T,Q,e]=solveQy(d,vs,ps,qz,qx,ind)
R=eye(3);
T=zeros(3,1);
Q=zeros(4,1);
d1=d(1);
d2=d(2);
d3=d(3);
d4=d(4);
d5=d(5);
d6=d(6);
d7=d(7);
d8=d(8);
d9=d(9);
d10=d(10);
d11=d(11);
d12=d(12);
d13=d(13);
d14=d(14);
d15=d(15);
d16=d(16);
d17=d(17);
d18=d(18);
d19=d(19);
d20=d(20);
d21=d(21);
d22=d(22);
d23=d(23);
d24=d(24);
d25=d(25);
d26=d(26);
d27=d(27);
d28=d(28);
d29=d(29);
d30=d(30);
d31=d(31);
d32=d(32);
d33=d(33);
d34=d(34);
d35=d(35);
d36=d(36);
d37=d(37);
d38=d(38);
d39=d(39);
d40=d(40);
d41=d(41);
d42=d(42);
d43=d(43);
d44=d(44);
d45=d(45);
d46=d(46);
dd =[ d25, d27 + d12*qx + d26*qz, d5*qx^2 + d13*qx*qz + d14*qx + d28*qz^2 + d29*qz + d30, d31, d2*qx^3 + d6*qx^2*qz + d7*qx^2 + d15*qx*qz^2 + d16*qx*qz + d17*qx + d32*qz^3 + d33*qz^2 + d34*qz + d36, d37 + d18*qx + d35*qz, d1*qx^4 + d3*qx^3*qz + d4*qx^3 + d8*qx^2*qz^2 + d9*qx^2*qz + d10*qx^2 + d19*qx*qz^3 + d20*qx*qz^2 + d21*qx*qz + d23*qx + d38*qz^4 + d39*qz^3 + d40*qz^2 + d42*qz + d44, d11*qx^2 + d22*qx*qz + d24*qx + d41*qz^2 + d43*qz + d45, d46];
d1=dd(1);
d2=dd(2);
d3=dd(3);
d4=dd(4);
d5=dd(5);
d6=dd(6);
d7=dd(7);
d8=dd(8);
d9=dd(9);
cy=[ 2*d1*d6 - d2*d4, d2*d6 - 2*d3*d4 + 4*d1*d8, 3*d2*d8 - 3*d4*d5, 2*d3*d8 - 4*d4*d7 - d5*d6, d5*d8 - 2*d6*d7];
yr=roots(cy);
%
e=-1;
for i=1:size(yr,1)
    if abs(real(yr(i)))<1E-8
        continue;
    end
    if abs(real(yr(i)))>1.1
        continue;
    end
    if abs(imag(yr(i)))>1E-8
        continue;
    end
    qy=real(yr(i));
    %[qx,qy,qz]
    [R_,T_,qz_,e_]=recover(vs,ps,qz,qx,qy,ind);
    if e<0 || e>e_
        R=R_;
        T=T_;
        e=e_;
        Q=qz_;
    end
end
end
function [R,T,Q,e]=recover(vs,ps,qz,qx,qy,ind)
global Rs;
global Ts;
global Got;
global A;
global B;
global v;
n=size(vs,2);

q=[1,1,1,1];
q(ind)=[qx,qy,qz];
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



