function out = Det(M)
n=size(M,1);
ps=perms(1:n);
N=size(ps,1);
out=0;
for i=1:N
    [i,N]
    index=ps(i,:);
    dis=numberofdisorder(index);
    temp=(-1)^dis;
    for j=1:n
       temp=temp*M(j,index(j));
    end
    out=out+temp;
end
end

function out=numberofdisorder(index)
out=0;
n=size(index,2);
for k=1:(n-1)
    for t=(k+1):n
        if index(k)>index(t)
            out=out+1;
        end
    end
end
end