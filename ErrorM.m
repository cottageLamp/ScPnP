function [ e ] = ErrorM( Rt,R,Tt,T )
e=[-1;-1];
%%
M=size(size(R),2);
if M>2
    for m=1:size(R,3)
        R_=reshape(R(:,:,m),3,3);
        T_=reshape(T(:,m),3,1);
        er=-1;
        for i=1:3
            ei=acos(Rt(:,i)'*R_(:,i))*180/pi;
            if er < ei
                er=ei;
            end
        end
        if e(1)<0 || e(1)>er
            e(1)=er;
            e(2)=norm(Tt-T_)/norm(Tt)*100;
        end
    end
    return;
end
%%
for i=1:3
    ei=Rt(:,i)'*R(:,i);
    ei=acos(ei)*180/pi;
    if e(1)<0 || e(1)> ei
        e(1)=ei;
    end
end
e(2)=norm(Tt-T)/norm(Tt)*100;
end

