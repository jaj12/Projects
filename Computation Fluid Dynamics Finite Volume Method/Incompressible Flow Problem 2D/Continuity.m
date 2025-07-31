function [acp,anp,asp,aep,awp,bp]=Continuity(dxc,dyc,dy,dx,Ux,Uy,Rho,mdotw,mdote,mdotn,mdots,acx,acy,n,m)

%Boundaries need to be taken into account.

for i=2:n-2
    for j=2:m-2
        Dcx(i,j)=dx(i)*dy(j)/acx(i,j);
        DEx(i,j)=dx(i)*dy(j)/acx(i+1,j);
        DNx(i,j)=dx(i)*dy(j)/acx(i,j+1);
        DWx(i,j)=dx(i)*dy(j)/acx(i-1,j);
        DSx(i,j)=dx(i)*dy(j)/acx(i,j-1);
        
        Dcy(i,j)=dx(i)*dy(j)/acy(i,j);
        DEy(i,j)=dx(i)*dy(j)/acy(i+1,j);
        DNy(i,j)=dx(i)*dy(j)/acy(i,j+1);
        DWy(i,j)=dx(i)*dy(j)/acy(i-1,j);
        DSy(i,j)=dx(i)*dy(j)/acy(i,j-1);
        
        Dex(i,j)=inter(Dcx(i,j),DEx(i,j),ge(i,j));
        Dnx(i,j)=inter(Dcx(i,j),DNx(i,j),gn(i,j));
        Dwx(i,j)=inter(Dcx(i,j),DWx(i,j),gw(i,j));
        Dsx(i,j)=inter(Dcx(i,j),DSx(i,j),gs(i,j));
        
        Dey(i,j)=inter(Dcy(i,j),DEy(i,j),ge(i,j));
        Dny(i,j)=inter(Dcy(i,j),DNy(i,j),gn(i,j));
        Dwy(i,j)=inter(Dcy(i,j),DWy(i,j),gw(i,j));
        Dsy(i,j)=inter(Dcy(i,j),DSy(i,j),gs(i,j));
        
        Dele(i,j)=Dex*dy(j);%simplified due to cartesian grid.
        Delw(i,j)=-Dwx*dy(j);
        Dels(i,j)=-Dsy*dx(i);
        Deln(i,j)=Dny*dx(i);
        
        awp(i,j)=-Rho(i,j)*Delw(i,j);
        aep(i,j)=-Rho(i,j)*Dele(i,j);
        anp(i,j)=-Rho(i,j)*Deln(i,j);
        asp(i,j)=-Rho(i,j)*Dels(i,j);
        acp(i,j)=-(awp(i,j)+aep(i,j)+anp(i,j)+asp(i,j));
        
        bp(i,j)=-(mdots+mdote+mdotn+mdotw);
    end
end