function [Tc]=TDMA(Tc,a,an,as,ae,aw,b,n,m)
an(:,m)=0;
ae(n,:)=0;
            for i=1:n-1
                for j=1:m-1
                    if j==1 && i==1 %South west CV
                        d(i,j)=b(i,j)-ae(i,j)*Tc(i+1,j);
                        P(i,j)=-an(i,j)/a(i,j);
                        Q(i,j)=d(i,j)/a(i,j);
                        Tc(i,j)=P(i,j)*Tc(i,j+1)+Q(i,j);

                    elseif j==1 && i~=1 && i~=n-1 % South CV
                        d(i,j)=b(i,j)-ae(i,j)*Tc(i+1,j) -aw(i,j)*Tc(i-1,j);
                        P(i,j)=-an(i,j)/a(i,j);
                        Q(i,j)=d(i,j)/a(i,j);

                    elseif j==1 && i==n-1 %South East CV
                        d(i,j)=b(i,j)-aw(i,j)*Tc(i-1,j);
                        P(i,j)=-an(i,j)/a(i,j);
                        Q(i,j)=d(i,j)/a(i,j);

                    elseif j~=1 && i==1 % West
                        d(i,j)=b(i,j)-ae(i,j)*Tc(i+1,j);
                        P(i,j)=-an(i,j)/(a(i,j)+as(i,j)*P(i,j-1));
                        Q(i,j)=(d(i,j)-as(i,j)*Q(i,j-1))/(a(i,j)+as(i,j)*P(i,j-1));

                    elseif j~=1 && i==n-1 %East
                        d(i,j)=b(i,j)-aw(i,j)*Tc(i-1,j);
                        P(i,j)=-an(i,j)/(a(i,j)+as(i,j)*P(i,j-1));
                        Q(i,j)=(d(i,j)-as(i,j)*Q(i,j-1))/(a(i,j)+as(i,j)*P(i,j-1));

                    else
                        d(i,j)=b(i,j)-aw(i,j)*Tc(i-1,j)-ae(i,j)*Tc(i+1,j);
                        P(i,j)=-an(i,j)/(a(i,j)+as(i,j)*P(i,j-1));
                        Q(i,j)=(d(i,j)-as(i,j)*Q(i,j-1))/(a(i,j)+as(i,j)*P(i,j-1));
                    end
                
                if j==m-1
                    Tc(i,j)=Q(i,j);
                else
                    Tc(i,j)=P(i,j)*Tc(i,j+1)+Q(i,j);
                end
                end
            end
end