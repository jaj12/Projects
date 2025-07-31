function [mdotn,mdotw,mdote,mdots]=RhieChow(ac,P,dxc,dyc,dx,dy,ge,gn,Rho,dpcxx,dpcyy,Ux,Uy)

n=length(dxc);
m=length(dyc);
    for i=1:n-2
        for j=1:m-2
        Dc(i,j)=dx(i)*dy(j)/ac(i,j);
        DE(i,j)=dx(i)*dy(j)/ac(i+1,j);
        DN(i,j)=dx(i)*dy(j)/ac(i,j+1);
        
        De(i,j)=inter(Dc(i,j),DE(i,j),ge(i,j));
        Dn(i,j)=inter(Dc(i,j),DN(i,j),gn(i,j));
        
        dpce(i,j)= dpcxx(i+1,j);
        dpcn(i,j)= dpcyy(i,j+1);
               
        
        dpe_(i,j)=inter(dpcxx(i,j),dpce(i,j),ge(i,j));
        dpn_(i,j)=inter(dpcyy(i,j),dpcn(i,j),gn(i,j));
        
        Ue_(i,j)=inter(Ux(i,j),Ux(i+1,j),ge(i,j));
        Un_(i,j)=inter(Uy(i,j),Uy(i,j+1),gn(i,j));
        
        dpe(i,j)=(P(i+1,j)-P(i,j))/dxc(i);
        dpn(i,j)=(P(i,j+1)-P(i,j))/dyc(j);
        
        Un=Un_(i,j)-Dn(i,j)*(dpn(i,j)-dpn_(i,j));
        Ue=Ue_(i,j)-De(i,j)*(dpe(i,j)-dpe_(i,j));
        end
    end
    
    for i=2:n-2
        for j=2:m-1
        Uw(i,j)=Ue(i-1,j);
        Us(i,j)=Un(i-1,j);
        mdotw(i,j)=-Uw(i,j)*Rho(i,j)*dy(j);
        mdotn(i,j)=Un(i,j)*Rho(i,j)*dx(i);
        mdote(i,j)=Ue(i,j)*Rho(i,j)*dy(j);
        mdots(i,j)=-Us(i,j)*Rho(i,j)*dx(i);
        end
    end
end

