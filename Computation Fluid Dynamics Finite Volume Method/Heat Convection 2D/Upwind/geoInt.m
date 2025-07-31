function [dx,dy,dxc,dyc,xc,yc,V,n,m]=geoInt(x,y)
n=length(x);
m=length(y);
dx(1:n-1)=x(2:n)-x(1:n-1);
dy(1:m-1)=y(2:m)-y(1:m-1);
%x-center and distance calculation

    for i=2:n
        for j=2:m

            V(i-1,j-1)=dx(i-1)*dy(j-1);
            if i~=2

            V(i-2,j-1)=dx(i-2)*dy(j-1);
            else
            if j==2
                yc(j-1)=dy(j-1)/2;
                dyc(j-1)=yc(j-1);
            else
                yc(j-1)=yc(j-2)+dy(j-1);
                dyc(j-1)=yc(j-1)-yc(j-2);
            end
            end
            
            if j~=2
                
            V(i-1,j-2)=dx(i-1)*dy(j-2);
            end

        end
        if i==2

            xc(i-1)=dx(i-1)/2;
            dxc(i-1)=xc(i-1);
        else

            xc(i-1)=xc(i-2)+dx(i-1);
            dxc(i-1)=xc(i-1)-xc(i-2);
        end
    end
    dyc(m)=y(end)-yc(m-1);
    dxc(n)=x(end)-xc(n-1);
end