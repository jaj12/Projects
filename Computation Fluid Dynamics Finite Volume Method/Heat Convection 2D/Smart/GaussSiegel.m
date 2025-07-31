function [Tc]=GaussSiegel(Tc,a,an,as,ae,aw,b,n,m)
for i=1:n-1
    for j=1:m-1

                    if i==1 && j==1 %South West Corner
                        Tc(i,j)=1/a(i,j)*(b(i,j)-an(i,j)*Tc(i,j+1)-ae(i,j)*Tc(i+1,j));

                    elseif i==1 && j==m-1 % North West Corner
                        Tc(i,j)=1/a(i,j)*(b(i,j)-as(i,j)*Tc(i,j-1)-ae(i,j)*Tc(i+1,j));  

                    elseif i==n-1 && j==1 % South East Corner
                        Tc(i,j)=1/a(i,j)*(b(i,j)-an(i,j)*Tc(i,j+1)-aw(i,j)*Tc(i-1,j));

                    elseif i==n-1 && j==m-1 % North East
                        Tc(i,j)=1/a(i,j)*(b(i,j)-as(i,j)*Tc(i,j-1)-aw(i,j)*Tc(i-1,j));

                    elseif i==n-1 && j~=1 && j~=m-1 % East Border
                        Tc(i,j)=1/a(i,j)*(b(i,j)-an(i,j)*Tc(i,j+1)-as(i,j)*Tc(i,j-1)-aw(i,j)*Tc(i-1,j));

                    elseif i==1 && j~=1 && j~=m-1 %West borner
                        Tc(i,j)=1/a(i,j)*(b(i,j)-an(i,j)*Tc(i,j+1)-as(i,j)*Tc(i,j-1)-ae(i,j)*Tc(i+1,j));

                    elseif j==1 && i~=1 && i~=n-1 % South Border
                        Tc(i,j)=1/a(i,j)*(b(i,j)-an(i,j)*Tc(i,j+1)-aw(i,j)*Tc(i-1,j)-ae(i,j)*Tc(i+1,j));

                    elseif j==m-1 && i~=1 && i~=n-1 % North Border
                        Tc(i,j)=1/a(i,j)*(b(i,j)-as(i,j)*Tc(i,j-1)-aw(i,j)*Tc(i-1,j)-ae(i,j)*Tc(i+1,j));

                    else
                        Tc(i,j)=1/a(i,j)*(b(i,j)-an(i,j)*Tc(i,j+1)-as(i,j)*Tc(i,j-1)-aw(i,j)*Tc(i-1,j)-ae(i,j)*Tc(i+1,j));
                        
                    end
    end
end
end
