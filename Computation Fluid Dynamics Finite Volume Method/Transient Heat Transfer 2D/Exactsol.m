function [T10]=Exactsol(Tmax,Tplot,fignum)
a=0.1;
b=0.2;
dt=1;
dx=0.0025;
dy=0.0025;
y=linspace(0.0025,b-0.0025,40);
x=linspace(0.0025,a-0.0025,20);
t=linspace(0,Tmax,Tmax+1);
T=zeros(length(x),length(y),length(t));
alpha=15/(500*7800);
m=0;
for k=1:length(t)
    for i=1:length(x)
        for j=1:length(y)
            for n=1:1000
                A0n=2000/(n*pi)*(1-cos(n*pi));
                T(i,j,k)=T(i,j,k)+A0n*cos(m*pi/a*x(i))*sin(n*pi/b*y(j))*exp(-alpha*pi^2*(m^2/a^2+n^2/b^2)*t(k));


            end
        end
    end
end


T10=T(:,:,Tplot+1);

figure(fignum)
namenum = num2str(Tplot);
contourf(x,y,T10')
 title(sprintf('Exact Temperature(K)t= %s s',namenum))
 xlabel('Length')
 ylabel('Width')
 color = colorbar; colormap(jet);
 ylabel(color, 'Temperature(K)')
end