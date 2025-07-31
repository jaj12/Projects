function [a,an,as,aw,ae,b]=Temperature(Tc,dxc,dyc,dy,dx,Qg,BE,BW,BN,BS,hbE,hbW,hbN,hbS,TinfE,TinfW,TinfN,TinfS,TbE,TbW,TbN,TbS,QbE,QbW,QbN,QbS,kn,ks,ke,kw,Ux,Uy,ConvectionScheme,Rho,TWin,TNin,TSin,TEin,Cp,dt,TransientScheme,Tcoldtrans,UNyin,USyin,UExin,UWxin,Uwx,Uex,Uny,Usy,n2i,s2i,e2i,w2i)
kbn= cond(TbN);
kbs= cond(TbS);
kbe= cond(TbE);
kbw= cond(TbW);

n=length(dxc);
m=length(dyc);
    for i=1:n-1
        for j=1:m-1
            
            %zeroing variables
            mdotw=0;
            mdots=0;
            mdotn=0;
            mdote=0;
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
            
            %changing boundary condition based on position.
            if i<n2i
                BNi=BN(1);
            else
                BNi=BN(2);
            end
            
            if i<s2i
                BSi=BS(1);
            else
                BSi=BS(2);
            end
            
            if j<e2i
                BEi=BE(1);
            else
                BEi=BE(2);
            end
            
            if j<w2i
                BWi=BW(1);
            else
                BWi=BW(2);
            end
            
            if i==1 %West
                
                if BWi==1
                ab=-kbw(j)*(dy(j)/(dxc(i)));
                bbw=-ab*TbW(j);
                awt=ab;

                elseif BWi==2
                bbw=QbW*dy(j);

                elseif BWi==3
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

                if BEi==1
                
                ab=-kbe(j)*(dy(j)/(dxc(i+1)));
                bbe=-ab*TbE(j);
                aet=ab;
                elseif BEi==2

                bbe=QbE*dy(j);

                elseif BEi==3
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
                
                if BSi==1

                ab=-kbs(i)*(dx(i)/(dyc(j)));
                bbs=-ab*TbS(i) + Qg*dx(i)*dy(j);
                ast=ab;
                elseif BSi==2

                bbs=QbS*dx(i);
                elseif BSi==3

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
                
                    if BNi==1

                    ab=-kbn(i)*(dx(i)/(dyc(j+1)));
                    bbn=-ab*TbN(i);
                    ant=ab;
                    elseif BNi==2
                    bbn=QbN*dx(i);

                    elseif BNi==3
                    RbN=(hbN*kbn(i)/dyc(j+1))/(hbN+kbn(i)/dyc(j+1));
                    ab=-RbN*dx(i);
                    bbn=-ab*TinfN;
                    ant=ab;
                    end
                    
            else
                   an(i,j)=-kn(i,j)*(dx(i)/dyc(j+1));
                   ant=an(i,j);
            end
          %%  Convection Coefficients

            bbwc=0;
            bbec=0;
            bbnc=0;
            bbsc=0;
            
            if i~=1
                mdotw=-Rho(i,j)*dy(j)*Uwx(i,j);
            else
                
            end
            
            if i~=n-1
                mdote=Rho(i,j)*dy(j)*Uex(i,j);
            end
            
            if j~=1
                mdots=-Rho(i,j)*dx(i)*Usy(i,j);
            end
            
            if j~=m-1
                mdotn=Rho(i,j)*dx(i)*Uny(i,j);
            end
           
            
            
                if i==1 %West
                    if BWi==4
                        mdotw=-Rho(i,j)*dy(j)*UWxin;
                        bbwc=-mdotw*Cp*TWin;
                        awt=awt+mdotw*Cp;
                    elseif BWi==5
                        mdotw=-Rho(i,j)*dy(j)*Ux(i,j);
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
                    else
                        Tcdimw=0;
                    end
                    Tfdimw=Smart(Tcdimw);
                    Tfw=Tfdimw*Qw+Tuw;
                    bdcw=mdotw*Cp*(Tfw-Tcw);
                end
                
                if j==1 %South
                    if BSi==4
                        mdots=-Rho(i,j)*dx(i)*USyin;
                        bbsc=-mdots*Cp*TSin;
                        ast=ast+Cp*mdots;
                    elseif BSi==5
                        bbsc=0;
                        mdots=-Rho(i,j)*dx(i)*Uy(i,j);
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
                    else
                        Tcdims=0;
                    end
                    Tfdims=Smart(Tcdims);
                    Tfs=Tfdims*Qs+Tus;
                    bdcs=mdots*Cp*(Tfs-Tcs);
                end
                
                if i==n-1 %East
                    if BEi==4
                        mdote=Rho(i,j)*dy(j)*UExin;
                        bbec=-Cp*mdote*TEin;
                        aet=aet+Cp*mdote;
                    elseif BEi==5
                        bbec=0;
                        mdote=Rho(i,j)*dy(j)*Ux(i,j);
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
                    else
                        Tcdime=0;
                    end
                    Tfdime=Smart(Tcdime);
                    Tfe=Tfdime*Qe+Tue;
                    bdce=mdote*Cp*(Tfe-Tce);
                end
                
                if j==m-1 %North
                    if BNi==4
                        mdotn=Rho(i,j)*dx(i)*UNyin;
                        bbnc=-Cp*mdotn*TNin;
                        ant=ant+Cp*mdotn;
                    elseif BNi==5
                        bbnc=0;
                        mdotn=Rho(i,j)*dx(i)*Uy(i,j);
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
                    else
                        Tcdimn=0;
                    end
                    Tfdimn=Smart(Tcdimn);
                    Tfn=Tfdimn*Qn+Tun;
                    bdcn=mdotn*Cp*(Tfn-Tcn);
                end
                
            bbctot=bbsc+bbnc+bbec+bbwc; %convection boundary
            bbtot=bbn+bbs+bbe+bbw; %diffusion boundary

                    
                
            if ConvectionScheme=="Upwind"
                bdctot=0;
            elseif ConvectionScheme=="Smart"
            bdctot=bdce+bdcw+bdcn+bdcs;%deffered correction
            end
            
            if TransientScheme=="BE" 
                a0=Rho(i,j)*Cp*dx(i)*dy(j)/dt;
            elseif TransientScheme=="CN"
                a0=Rho(i,j)*Cp*dx(i)*dy(j)/(dt/2);
            else
                a0=0;
            end

                
            a(i,j)=-(ant+ast+aet+awt)+Cp*(mdote+mdotw+mdots+mdotn)+a0;
            b(i,j)=-bdctot+bbctot+bbtot+Qg*dx(i)*dy(j)+a0*Tcoldtrans(i,j);
        end
    end
end