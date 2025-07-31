function [gn,gs,ge,gw]= Geo(V)

[x,y]=size(V);
    for i=1:x
        for j=1:y
            if j==1
                gn(i,j)= V(i,j)/(V(i,j)+V(i,j+1));
            elseif j==y
                gs(i,j)= V(i,j)/(V(i,j)+V(i,j-1));
            else
                gn(i,j)= V(i,j)/(V(i,j)+V(i,j+1));
                gs(i,j)= V(i,j)/(V(i,j)+V(i,j-1));
            end
            if i==1
                ge(i,j)= V(i,j)/(V(i,j)+V(i+1,j));
                
            elseif i==x
                gw(i,j)= V(i,j)/(V(i,j)+V(i-1,j));
                
            else
                gw(i,j)= V(i,j)/(V(i,j)+V(i-1,j));
                ge(i,j)= V(i,j)/(V(i,j)+V(i+1,j));
            end
        end
    end
end