function Rscaled=ResidCalc(a,aw,an,as,ae,Tc,b,n,m)
an(:,m)=0;
ae(n,:)=0;

sum=0;
for i=1:n-1
    for j=1:m-1
        
        sum(i,j)=sum(i,j)+a(i,j)*Tc(i,j)-b(i,j);
        if i~=1
        sum(i,j)=sum(i,j)+aw(i,j)*Tc(i-1,j);
        end
        if j~=1
        sum(i,j)=sum(i,j)+as(i,j)*Tc(i,j-1);
        end
        if i~=n-1
        sum(i,j)=sum(i,j)+ae(i,j)*Tc(i+1,j);
        end
        if j~=m-1
        sum(i,j)=sum(i,j)+an(i,j)*Tc(i,j+1);
        end

    end
end

scale=max(ac.*Tc);

Rscaled=max(sum/scale);

end