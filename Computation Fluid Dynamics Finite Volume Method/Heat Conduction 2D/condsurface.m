function [kn,ks,ke,kw]=condsurface(kc,gn,gs,ge,gw)

x=length(gn);
y=size(ge,2);
    for i=1:x
        for j=1:y
            if j==1
                kn(i,j)= inter(kc(i,j),kc(i,j+1),gn(i,j));
                
            elseif j==y
                ks(i,j)= inter(kc(i,j),kc(i,j-1),gs(i,j));
            else
                kn(i,j)= inter(kc(i,j),kc(i,j+1),gn(i,j));
                ks(i,j)= inter(kc(i,j),kc(i,j-1),gs(i,j));
            end
            if i==1
                ke(i,j)= inter(kc(i,j),kc(i+1,j),ge(i,j));
            elseif i==x
                
                kw(i,j)= inter(kc(i,j),kc(i-1,j),gw(i,j));
            else
                
                kw(i,j)= inter(kc(i,j),kc(i-1,j),gw(i,j));
                ke(i,j)= inter(kc(i,j),kc(i+1,j),ge(i,j));
            end
        end
    end
end
