function [Pn,Ps,Pe,Pw]=Pressuresurface(Pc,gn,gs,ge,gw)

x=size(gn,1);
y=size(ge,2);
    for i=1:x
        for j=1:y
            if j==1
                Pn(i,j)= inter(Pc(i,j),Pc(i,j+1),gn(i,j));
                
            elseif j==y
                Ps(i,j)= inter(Pc(i,j),Pc(i,j-1),gs(i,j));
            else
                Pn(i,j)= inter(Pc(i,j),Pc(i,j+1),gn(i,j));
                Ps(i,j)= inter(Pc(i,j),Pc(i,j-1),gs(i,j));
            end
            if i==1
                Pe(i,j)= inter(Pc(i,j),Pc(i+1,j),ge(i,j));
            elseif i==x
                
                Pw(i,j)= inter(Pc(i,j),Pc(i-1,j),gw(i,j));
            else
                Pw(i,j)= inter(Pc(i,j),Pc(i-1,j),gw(i,j));
                Pe(i,j)= inter(Pc(i,j),Pc(i+1,j),ge(i,j));
            end
        end
    end
end