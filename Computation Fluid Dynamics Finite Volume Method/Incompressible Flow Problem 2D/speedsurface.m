function [Uwx,Uex,Uny,Usy]=speedsurface(Ux,Uy,gn,gs,ge,gw)

x=size(gn,1);
y=size(ge,2);
    for i=1:x
        for j=1:y
            if j==1
                Uny(i,j)= inter(Uy(i,j),Uy(i,j+1),gn(i,j));
                
            elseif j==y
                Usy(i,j)= inter(Uy(i,j),Uy(i,j-1),gs(i,j));
            else
                Uny(i,j)= inter(Uy(i,j),Uy(i,j+1),gn(i,j));
                Usy(i,j)= inter(Uy(i,j),Uy(i,j-1),gs(i,j));
            end
            if i==1
                Uex(i,j)= inter(Ux(i,j),Ux(i+1,j),ge(i,j));
            elseif i==x
                
                Uwx(i,j)= inter(Ux(i,j),Ux(i-1,j),gw(i,j));
            else
                Uwx(i,j)= inter(Ux(i,j),Ux(i-1,j),gw(i,j));
                Uex(i,j)= inter(Ux(i,j),Ux(i+1,j),ge(i,j));
            end
        end
    end
end