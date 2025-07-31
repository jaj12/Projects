function [a,an,as,aw,ae,b]=CoeffCal(Tc,dxc,dyc,dy,dx,Qg,BE,BW,BN,BS,hbE,hbW,hbN,hbS,TinfE,TinfW,TinfN,TinfS,TbE,TbW,TbN,TbS,QbE,QbW,QbN,QbS,kn,ks,ke,kw,Ux,Uy,ConvectionScheme,Rho,TWin,TNin,TSin,TEin,Cp)
kbn= cond(TbN);
kbs= cond(TbS);
kbe= cond(TbE);
kbw= cond(TbW);

n=length(dxc);
m=length(dyc);
    for i=1:n-1
        for j=1:m-1
        
            bdce=0;
            bdcw=0;
            bdcn=0;
            bdcs=0;
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

            bbwc=0;
            bbec=0;
            bbnc=0;
            bbsc=0;
            
            mdote=Rho(i,j)*dy(j)*Ux(i,j);
            mdotw=-Rho(i,j)*dy(j)*Ux(i,j);
            mdots=-Rho(i,j)*dx(i)*Uy(i,j);
            mdotn=Rho(i,j)*dx(i)*Uy(i,j);
            
            
                if i==1 %West
                    if BW==4
                        bbwc=-mdotw*Cp*TWin;
                        awt=awt+mdotw*Cp;
                    elseif BW==5
                        bbwc=0;
                    end
                else
                    if Ux>0
                    if i==2
                        Tuw=TWin;
                    else
                        Tuw=Tc(i-2,j);
                    end
                    Tdw=Tc(i,j);
                    Tcw=Tc(i-1,j);
                    awc=Cp*mdotw;
                    aw(i,j)=aw(i,j)+awc;
                    awt=aw(i,j);
                    else
                        if i==n-1
                            Tuw=TEin;
                        else
                            Tuw=Tc(i+1,j);
                        end
                        Tdw=Tc(i-1,j);
                        Tcw=Tc(i,j);
                    end
                    Qw=Tdw-Tuw;
                    if Qw~=0
                        Tcdimw=(Tcw-Tuw)/Qw;
                        Tfdimw=Smart(Tcdimw);
                        Tfw=Tfdimw*Qw+Tuw;
                    else
                        Tfw=Tcw;
                    end
                    
                    bdcw=mdotw*Cp*(Tfw-Tcw);
                end
                
                if j==1 %South
                    if BS==4
                        bbsc=-mdots*Cp*TSin;
                        ast=ast+Cp*mdots;
                    elseif BS==5
                        bbsc=0;
                    end
                else
                    if Uy>0
                    if j==2
                        Tus=TSin;
                    else
                        Tus=Tc(i,j-2);
                    end
                    Tds=Tc(i,j);
                    Tcs=Tc(i,j-1);
                    asc=Cp*mdots;
                    as(i,j)=as(i,j)+asc;
                    ast=as(i,j);
                    else
                        if j==m-1
                            Tus=TNin;
                        else
                            Tus=Tc(i,j+1);
                        end
                        
                        Tds=Tc(i,j-1);
                        Tcs=Tc(i,j);
                    end
                    Qs=Tds-Tus;
                    if Qs~=0
                        Tcdims=(Tcs-Tus)/Qs;
                        Tfdims=Smart(Tcdims);
                        Tfs=Tfdims*Qs+Tus;
                    else
                        Tfs=Tcs;
                    end
                    
                    bdcs=mdots*Cp*(Tfs-Tcs);
                end
                
                if i==n-1 %East
                    if BE==4
                        bbec=-Cp*mdote*TEin;
                        aet=aet+Cp*mdote;
                    elseif BE==5
                        bbec=0;
                    end
                else
                    if Ux<0
                        aec=Cp*mdote;
                        ae(i,j)=ae(i,j)+aec;
                        aet=ae(i,j);
                        if i==n-2
                            Tue=TEin;
                        else
                            Tue=Tc(i+2,j);
                        end
                        Tde=Tc(i,j);
                        Tce=Tc(i+1,j);
                    else
                      
                        if i==1
                            Tue=TWin;
                        else
                            Tue=Tc(i-1,j);
                        end
                        Tde=Tc(i+1,j);
                        Tce=Tc(i,j);
                    end
                    Qe=Tde-Tue;
                    if Qe~=0
                        Tcdime=(Tce-Tue)/Qe;
                        Tfdime=Smart(Tcdime);
                        Tfe=Tfdime*Qe+Tue;
                    else
                        Tfe=Tce;
                    end
                    
                    bdce=mdote*Cp*(Tfe-Tce);
                end
                
                if j==m-1 %North
                    if BN==4
                        bbnc=-Cp*mdotn*TNin;
                        ant=ant+Cp*mdotn;
                    elseif BN==5
                        bbnc=0;
                    end
                else
                    if Uy<0
                        anc=Cp*mdotn;
                        an(i,j)=an(i,j)+anc;
                        ant=an(i,j);
                        if n==n-2
                            Tun=TNin;
                        else
                            Tun=Tc(i,j+1);
                        end
                        Tdn=Tc(i,j);
                        Tcn=Tc(i,j+1);
                    else
                        if j==1
                            Tun=TSin;
                        else
                            Tun=Tc(i,j-1);
                        end
                        Tdn=Tc(i,j+1);
                        Tcn=Tc(i,j);
                    end
                    Qn=Tdn-Tun;
                    if Qn~=0
                        Tcdimn=(Tcn-Tun)/Qn;
                        Tfdimn=Smart(Tcdimn);
                        Tfn=Tfdimn*Qn+Tun;
                    else
                        Tfn=Tcn;
                    end
                    
                    bdcn=mdotn*Cp*(Tfn-Tcn);
                end
                
            bbctot=bbsc+bbnc+bbec+bbwc; %convection boundary
            bbtot=bbn+bbs+bbe+bbw; %diffusion boundary

            
            
            if ConvectionScheme=="Upwind"
                bdctot=0;
            elseif ConvectionScheme=="Smart"
            bdctot=bdce+bdcw+bdcn+bdcs;%deffered correction
            end

            a(i,j)=-(ant+ast+aet+awt)+Cp*(mdote+mdotw+mdots+mdotn);
            b(i,j)=-bdctot+bbctot+bbtot+Qg*dx(i)*dy(j);
        end
    end
end