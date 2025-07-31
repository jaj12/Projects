function [a,an,as,aw,ae,b]=CoeffCal(Tc,dxc,dyc,dy,dx,Qg,BE,BW,BN,BS,hbE,hbW,hbN,hbS,TinfE,TinfW,TinfN,TinfS,TbE,TbW,TbN,TbS,QbE,QbW,QbN,QbS,kn,ks,ke,kw,Ux,Uy,ConvectionScheme,Rho,UxinW,UxinE,UyinS,UyinN,TWin,TNin,TSin,TEin,Cp)
kbn= cond(TbN);
kbs= cond(TbS);
kbe= cond(TbE);
kbw= cond(TbW);

n=length(dxc);
m=length(dyc);
    for i=1:n-1
        for j=1:m-1
        
            
            ant=0;
            ast=0;
            awt=0;
            aet=0;
            bbw=0;
            bbs=0;
            bbn=0;
            bbe=0;
            if i==1 %West
                
                if BW==1
                ab=-kbw(j)*(dy(j)/(dxc(i)));
                bbw=-ab*TbW(j);
                awt=ab;

                elseif BW==2
                bbw=QbW*dy(j);

                elseif BW==3
                RbW=(hbW*kbw(j)/dxc(i))/(hbW+kbw/dxc(i));
                ab=-RbW*dy(j);
                awt=ab;
                bbw=-ab*TinfW;
                end
                
            else
                aw(i,j)=-kw(i,j)*(dy(j)/dxc(i));
                awt=aw(i,j);
            end
            
            if i==n-1 %East

                if BE==1
                
                ab=-kbe(j)*(dy(j)/(dxc(i+1)));
                bbe=-ab*TbE(j);
                aet=ab;
                elseif BE==2

                bbe=QbE*dy(j);

                elseif BE==3
                RbE=(hbE*kbe(j)/dxc(i+1))/(hbE+kbe(j)/dxc(i+1));
                ab=-RbE*dy(j);
                bbe=-ab*TinfE;
                aet=ab;
                end
                
            else
                ae(i,j)=-ke(i,j)*(dy(j)/dxc(i+1));
                aet=ae(i,j);
            end
                
            if j==1 % South
                
                if BS==1

                ab=-kbs(i)*(dx(i)/(dyc(j)));
                bbs=-ab*TbS(i) + Qg*dx(i)*dy(j);
                ast=ab;
                elseif BS==2

                bbs=QbS*dx(i);
                elseif BS==3

                RbS=(hbS*kbs(i)/dyc(j))/(hbS+kbs(i)/dyc(j));
                ab=-RbS*dx(i);
                bbs=-ab*TinfS;
                ast=ab;
                end
                
            else
                as(i,j)=-ks(i,j)*(dx(i)/dyc(j));
                ast=as(i,j);
            end
            if j==m-1 % North
                    if BN==1

                    ab=-kbn(i)*(dx(i)/(dyc(j+1)));
                    bbn=-ab*TbN(i);
                    ant=ab;
                    elseif BN==2
                    bbn=QbN*dx(i);

                    elseif BN==3
                    RbN=(hbN*kbn(i)/dyc(j+1))/(hbN+kbn(i)/dyc(j+1));
                    ab=-RbN*dx(i);
                    bbn=-ab*TinfN;
                    ant=ab;
                    end
                    
            else
                   an(i,j)=-kn(i,j)*(dx(i)/dyc(j+1));
                   ant=an(i,j);
            end
            %%Convection Coefficients
            
            aec=0;
            awc=0;
            anc=0;
            asc=0;
            bbwc=0;
            bbec=0;
            bbnc=0;
            bbsc=0;
            if ConvectionScheme=="Upwind"
            
                if i==1 %West
                    if BW==4
                        bbwc=Cp*Rho(i,j)*TWin*dy(j)*UxinW;
                        awt=-Cp*Rho(i,j)*dy(j)*UxinW;
                    elseif BW==5
                        bbwc=0;
                    end
                    
                else
                    
                    if Ux>0
                        awc=-Cp*Rho(i,j)*Ux(i,j)*dy(j);
                        aw(i,j)=aw(i,j)+awc;
                        awt=aw(i,j);
                    end
                end
                
                if j==1 %South
                    if BS==4
                        bbsc=Cp*Rho(i,j)*TSin*dx(i)*UyinS;
                        ast=-Cp*Rho(i,j)*dx(i)*UyinS;
                    elseif BS==5
                        bbsc=0;
                    end
                    
                else
                    
                    if Uy>0
                        asc=-Cp*Rho(i,j)*Uy(i,j)*dx(i);
                        as(i,j)=as(i,j)+asc;
                        ast=as(i,j);
                    end
                end
                
                if i==n-1 %East
                    if BE==4
                        bbec=Cp*Rho(i,j)*TEin*(dy(j)*UxinE);
                        aet=-Cp*Rho(i,j)*(dy(j)*UxinE);
                    elseif BE==5
                        bbec=0;
                    end
                    
                else
                    
                    if Ux<0
                        aec=-Cp*Rho(i,j)*(-Ux(i,j))*dy(j);
                        ae(i,j)=ae(i,j)+aec;
                        aet=ae(i,j);
                    end
                
                end
                
                if j==m-1 %North
                    if BN==4
                        bbnc=Cp*Rho(i,j)*TNin*(dx(i)*UyinN);
                        ant=-Cp*Rho(i,j)*dx(i)*UyinN;
                    elseif BN==5
                        bbnc=0;
                    end
                    
                else
                    
                    if Uy<0
                        anc=-Cp*Rho(i,j)*(-Uy(i,j))*dx(i);
                        an(i,j)=an(i,j)+anc;
                        ant=an(i,j);
                    end
                end
                
            end
            
            
            
            bconvection=bbnc+bbsc+bbec+bbwc;
            bdiffusion=bbn+bbs+bbe+bbw;
            
            a(i,j)=-(ant+ast+aet+awt);
            b(i,j)=bdiffusion+bconvection+Qg*dx(i)*dy(j);
        end
    end
end