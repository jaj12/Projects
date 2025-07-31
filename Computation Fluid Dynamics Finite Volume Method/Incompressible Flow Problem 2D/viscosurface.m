function [mun,mus,mue,muw]=viscosurface(muc,gn,gs,ge,gw)

x=size(gn,1);
y=size(ge,2);
    for i=1:x
        for j=1:y
            if j==1
                mun(i,j)= inter(muc(i,j),muc(i,j+1),gn(i,j));
                
            elseif j==y
                mus(i,j)= inter(muc(i,j),muc(i,j-1),gs(i,j));
            else
                mun(i,j)= inter(muc(i,j),muc(i,j+1),gn(i,j));
                mus(i,j)= inter(muc(i,j),muc(i,j-1),gs(i,j));
            end
            if i==1
                mue(i,j)= inter(muc(i,j),muc(i+1,j),ge(i,j));
            elseif i==x
                
                muw(i,j)= inter(muc(i,j),muc(i-1,j),gw(i,j));
            else
                muw(i,j)= inter(muc(i,j),muc(i-1,j),gw(i,j));
                mue(i,j)= inter(muc(i,j),muc(i+1,j),ge(i,j));
            end
        end
    end
end