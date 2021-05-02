reset(symengine);
clc;
clear;

% the unknowns
syms qw qx qy qz x4;
R=[qw^2+qx^2 - qy^2 - qz^2,2*qx*qy - 2*qz*qw,2*qx*qz + 2*qy*qw;...
2*qx*qy + 2*qz*qw,qw^2 + qy^2 - qx^2 - qz^2,2*qy*qz - 2*qx*qw;...
2*qx*qz - 2*qy*qw,2*qy*qz + 2*qx*qw,qw^2 + qz^2 - qx^2 - qy^2];

% a point that is given by measurements
syms vx vy vz px py pz;
vi=[vx;vy;vz];
viT=[vx,vy,vz];
pi=[px;py;pz];
piT=[px,py,pz];

% intermedia variable
syms v1 v2 v3;
v=[v1;v2;v3]; % v=\sum vi
vT=[v1,v2,v3];

S_ViRpi=vi*viT*R*pi; % \sum (I-vi*vi')R*pi = \sum -vi*vi'*R*pi
[c_,t_]=coeffs(S_ViRpi(1),[qw,qx,qy,qz]);
syms a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11; 
% expressions of a1, ..., a11 are given in c_
% for example a1 =\sum c_(1)
S_ViRpi(1)=sum([a1 a2 a3 a4 a5 a6 a7 a8 a9 a10].*t_);
[c_,t_]=coeffs(S_ViRpi(2),[qw,qx,qy,qz]);
syms b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11; % similar to a1, ...., a11
S_ViRpi(2)=sum([b1 b2 b3 b4 b5 b6 b7 b8 b9 b10].*t_);
[c_,t_]=coeffs(S_ViRpi(3),[qw,qx,qy,qz]);
syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11;%  similar to a1, ...., a11
S_ViRpi(3)=sum([c1 c2 c3 c4 c5 c6 c7 c8 c9 c10].*t_);

S_viRpi=viT*R*pi;
[c_,t_]=coeffs(S_viRpi,[qw,qx,qy,qz]);
syms d1 d2 d3 d4 d5 d6 d7 d8 d9 d10 d11; % similar to a1, ...., a11
S_viRpi=sum([d1 d2 d3 d4 d5 d6 d7 d8 d9 d10].*t_);

% intermedia matrices. See our paper for their definition and expression
syms U11 U12 U13 U21 U22 U23 U31 U32 U33;
U=[U11,U12,U13;U21 U22 U23; U31 U32 U33];
syms W11 W12 W13 W21 W22 W23 W31 W32 W33;
W=[W11,W12,W13;W21 W22 W23; W31 W32 W33];

% U0 and U1 are temporary variables
U0=U*(W*S_ViRpi*x4+v*S_viRpi*x4-v);
[c_,t_]=coeffs(U0(1),[qw,qx,qy,qz,x4]);
syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11;% expressions of u1, ..., u11 are given in c_
U0(1)=sum([u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11].*t_);
[c_,t_]=coeffs(U0(2),[qw,qx,qy,qz,x4]);
syms v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11;% similar to u1, ..., u11
U0(2)=sum([v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11].*t_);
[c_,t_]=coeffs(U0(3),[qw,qx,qy,qz,x4]);
syms w1 w2 w3 w4 w5 w6 w7 w8 w9 w10 w11;% similar to u1, ..., u11
U0(3)=sum([w1 w2 w3 w4 w5 w6 w7 w8 w9 w10 w11].*t_);

U1=R*pi*x4+U0;

obj=sum(U1.*((eye(3)-vi*viT)*U1)); % note that (eye(3)-vi*vi')'*(eye(3)-vi*vi') =(eye(3)-vi*vi')
obj=subs(obj,qw,1); % may switch to qx qy or qz
[c_,t_]=coeffs(obj,[qx,qy,qz,x4]);
syms o1 o2 o3 o4 o5 o6 o7 o8 o9 o10 o11 o12 o13 o14 o15 o16 o17 o18 o19 o20 o21 o22 o23 o24 o25 o26 o27 o28 o29 o30 o31 o32 o33 o34 o35 o36 o37 o38 o39 o40 o41 o42 o43 o44 o45 o46;
% the expressions of o1, ..., o46 are given in c_
% for i=1:size(c_,2)
%    fprintf('%s,...\r',char(c_(i))); 
% end
os=[o1 o2 o3 o4 o5 o6 o7 o8 o9 o10 o11 o12 o13 o14 o15 o16 o17 o18 o19 o20 o21 o22 o23 o24 o25 o26 o27 o28 o29 o30 o31 o32 o33 o34 o35 o36 o37 o38 o39 o40 o41 o42 o43 o44 o45 o46];
obj=sum(os.*t_);
% obj is our objective function

% use 'diff' function to compute the optimal conditions
c1=diff(obj,qx);
[c1_,t1_]=coeffs(c1,x4);% t1_=[x4^2,x4] and it is equivalent to use [x4,1]
c1=c1_(1)*x4+c1_(2);

c2=diff(obj,qy);
[c2_,t2_]=coeffs(c2,x4);
c2=c2_(1)*x4+c2_(2);

c3=diff(obj,qz);
[c3_,t3_]=coeffs(c3,x4);
c3=c3_(1)*x4+c3_(2);

c4=diff(obj,x4);
[c4_,t4_]=coeffs(c4,x4);

% call our dixon function to compute the Dixion matrix
[M,tx]=dixon([c1,c2,c3,c4],[qx,qy,x4]);% solve for qz. One may switch to qx or qy

% check that our problem satisfies a condition given in 
% "Algebraic and geometric reasoning using dixon resultants, in: Proceedings of the international symposium on Symbolic and algebraic computation"
% M1=M;
% for i=1:size(M1,2)
%     index=setdiff(1:size(M1,2),i);
%     rank(M1(:,index))% exist < 22
% end

rows=setdiff(1:size(M,1),[1,2,3]);
cols=setdiff(1:size(M,2),[1,2,4,5,6,8,11,14,19]);
D=M(rows,cols); % remove linear dependent rows and get our Dixion matrix

% output expressions of the matrix polynomial of the Eq.(27) of our paper
M2NFile(D,'zMqz',qz);