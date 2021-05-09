function [ M,tx ] = dixon(fs,x)
n=size(fs,2);
m=size(x,2);

syms t1 t2 t3 t4;
t=[t1,t2,t3,t4];
t=t(1:m);

D=sym(zeros(n,n));
for i=1:n
   D(1,i)=fs(i); 
end
for i=2:n
   D(i,:)=subs(D(1,:),x(1:(i-1)),t(1:(i-1))); 
end
delta=Det(D);
ds=1;

for i=1:m
    ds=ds*(x(i)-t(i));
end

s=feval(symengine,'divide',delta,ds);
if s(2)~=0
    fprintf('dixon dividing error\n');
    return;
else
   delta=s(1);
end
fprintf('delta\n');

[~,tx]=coeffs(delta,x);
allTerms=sum(tx);
fprintf('all terms:%s\n',char(tx));

[C,~]=coeffs(delta,t);
fprintf('all functions\n');

syms a1;
M=sym(zeros(size(C,2),size(tx,2)));
fprintf('size of M: %dx%d\n',size(M));
for i=1:size(C,2)
    [cx,tx]=coeffs(C(i)+a1*allTerms,x);
    M(i,:)=cx-a1;
end
fprintf('rank of M: %d\n',rank(M));
end

