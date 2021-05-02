function  M2NFile( M ,name,xname)
m=size(M,2);

temp=0;
for i=1:m
    for j=1:m
        [c,t]=coeffs(M(i,j),xname);
        temp=temp+sum(t);
    end
end
[~,at]=coeffs(temp,xname);
n=size(at,2);
syms aa;
at=aa*sum(at);
% 
Ms=sym(zeros(n,m,m));
for i=1:m
    for j=1:m
        [c,t]=coeffs(M(i,j)+at,xname);
        for k=1:n
           Ms(k,i,j)=c(k)-aa; 
        end
    end
end
% An=reshape(Ms(n,:,:),[m,m]);
% rank(An)

for k=1:n
    filename=sprintf('%s%d.m',name,k);
    fid=fopen(filename,'w');
    fprintf(fid,'function out=%s%d(os)\r\n',name,k);
    for i=1:46
        fprintf(fid,'o%d=os(%d);\r',i,i,i,i,i,i);
    end
    fprintf(fid,'out=zeros(%d,%d);\r\n',m,m);
    for i=1:m
        for j=1:m
               %temp=simplify(Ms(k,i,j));
               temp=(Ms(k,i,j));
               if temp==0
                   continue;
               end
               temp=horner(temp);%Shorten(temp,cds);
               [i,j,k]
               fprintf(fid,'out(%d,%d)=%s;\r\n',i,j,char(temp)); 
            end
        end
    fprintf(fid,'end');
    fclose(fid);
end
end

