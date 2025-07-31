function [aw,an,ae,as,a,b,dpcx,dpcy,Pwb,Peb,Pnb,Psb]=VelocityY(P,dxc,dyc,dy,dx,BE,BW,BN,BS,WN,WE,WS,WW,Ux,Uy,ConvectionScheme,Rho,dt,TransientScheme,UNiny,USiny,UWiny,UEiny,UEinx,UWinx,Uwx,Uex,Uny,Usy,Uyoldtrans,TbN,TbS,TbE,TbW,n2i,s2i,e2i,w2i,PbW,PbE,PbS,PbN)

mbun= visco(TbN);
mbus= visco(TbS);
mbue= visco(TbE);
mbuw= visco(TbW);

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
            bbsy=0;
            bbny=0;
            bbex=0;
            bbwx=0;
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
            
            if i<n2i
                WNi=WN(1);
            else
                WNi=WN(2);
            end
            
            if i<s2i
                WSi=WS(1);
            else
                WSi=WS(2);
            end
            
            if j<e2i
                WEi=BW(1);
            else
                WEi=WE(2);
            end
            
            if j<w2i
                WWi=WW(1);
            else
                WWi=WW(2);
            end
            
            if i==1 %West
                dpcx(i,j)=(P(i+1,j)-PbW(j))/(dxc(i)+dxc(i+1));
                dUyw(i,j)=(Uy(i,j)-UbWy(j))/dxc(i);
                if WWi==9%no slip
                aby=-mbuw(j)*(dy(j)/(dxc(i)));
                Pwb(i)=P(i,j)-dpcx(i,j)*dxc(i);
                Pwt=Pwb(j);
                bbwx=Pwb(i)*dy(j)-aby*UWwally;
                
                awt=aby;
                elseif WWi==10 %slip
                Pwb(j)=P(i,j)-dpcx(i,j)*dxc(i);
                Pwt=Pwb(j);
                bbwx=Pwb(j)*dy(j);

                elseif BWi==4 %Inlet Specified velocity inlet
                Pwb(j)=P(i,j)-dpcx(i,j)*dxc(i);
                Pwt=Pwb(j);
                mdotw=-Rho(i,j)*dy(j)*UWinx;
                bbwx=-mdotw*UWiny+mbuw(j)*(dy(j)/(dxc(i)))*UWiny+Pwb(j)*dy(j);
                awt=awt+mdotw;
                elseif BWi==5 %Inlet Specified Static Pressure and Velocity Direction
                
                elseif BWi==6 %Inlet Specified Total Pressure and Velocity Direction
                    
                elseif BWi==7 %Outlet Specified Static Pressure
                mdotw=-Rho(i,j)*dy(j)*Ux(i,j);
                Pwb(j)=PWin;
                Pwt=PWin;
                bbwx=PWout*dy(j);
                elseif BWi==8 %Outlet Fully Developed Outlet Flow
                Pwb(j)=P(i,j)-dpcx(i,j)*dxc(i);
                Pwt=Pwb(j);
                mdotw=-Rho(i,j)*dy(j)*Ux(i,j);
                bbwx=-mdotw*dux(i,j)*dx(i+1)+Pwb(j)*dy(j);
                end
                
            else
                
                dpcx(i,j)=(P(i+1,j)-P(i-1,j))/(dxc(i)+dxc(i+1));
                dpcy(i,j)=(P(i,j+1)-P(i,j-1))/(dyc(j)+dyc(j+1));
                dUyw(i,j)=(Uy(i,j)-Uy(i-1,j))/dxc(i);
                aw(i,j)=-muw(i,j)*(dy(j)/dxc(i));
                awt=aw(i,j);
                Pnt=Pn(i,j);
                Pet=Pe(i,j);
                Pst=Ps(i,j);
            end
            
            if i==n-1 %East
                dpcx(i,j)=-(P(i-1,j)-PbE(j))/(dxc(i)+dxc(i+1));
                dUye(i,j)=-(Uy(i,j)-UbEy(j))/dxc(i+1);
                if WEi==9
                aby=-mbue(j)*(dy(j)/(dxc(i+1)));
                Peb(j)=P(i,j)+dpcx(i,j)*dxc(i+1);
                Pet=Peb(j);
                bbex=Peb(j)*dy(j)-aby*UEwally;
                
                aet=aby;
                elseif WEi==10
                Peb(j)=P(i,j)+dpcx(i,j)*dxc(i+1);
                Pet=Peb(j);
                bbex=Peb(j)*dy(j);
                

                elseif BEi==4
                Peb(j)=P(i,j)+dpcx(i,j)*dxc(i+1);
                Pet=Peb(j);
                mdote=Rho(i,j)*dy(j)*UEinx;
                bbex=-mdote*UEiny+mbue(j)*(dy(j)/(dxc(i+1)))*UEiny+Peb(j)*dy(j);
                aet=aet+mdote;
                elseif BEi==5
                    
                elseif BEi==6
                elseif BEi==7
                mdote=Rho(i,j)*dy(j)*Ux(i,j);
                Peb(j)=PEin;
                Pet=PEin;
                bbex=PEout*dy(j);
                elseif BEi==8
                Peb(j)=P(i,j)+dpcx(i,j)*dxc(i+1);
                Pet=Peb(j);
                mdote=Rho(i,j)*dy(j)*Ux(i,j);
                bbex=-mdote*dux(i,j)*dxc(i+1)+Peb(j)*dy(j);
                end
                
            else
                
                dpcx(i,j)=(P(i+1,j)-P(i-1,j))/(dxc(i)+dxc(i+1));
                dpcy(i,j)=(P(i,j+1)-P(i,j-1))/(dyc(j)+dyc(j+1));
                dUye(i,j)=(Uy(i+1,j)-Uy(i,j))/dxc(i+1);
                ae(i,j)=-mue(i,j)*(dy(j)/dxc(i+1));
                aet=ae(i,j);
                Pnt=Pn(i,j);
                Pwt=Pw(i,j);
                Pst=Ps(i,j);
            end
                
            if j==1 % South
                dpcy(i,j)=(P(i,j+1)-PbS(i))/(dyc(j)+dyc(j+1));
                dUys(i,j)=(Uy(i,j)-UbSy(i))/dyc(j);
                if WSi==9
                abx=-mbus(j)*(dx(i)/(dyc(j)));
                Psb(i)=P(i,j)-dpcy(i,j)*dyc(j);
                Pst=Psb(i);
                bbsy=Psb(j)*dx(i)-abx*USwally;
                
                ast=abx;
                elseif WSi==10
                Psb(i)=P(i,j)-dpcy*dyc(j);
                Pst=Psb(i);
                bbsy=Psb(j)*dx(i);
                
                elseif BSi==4
                Psb(i)=P(i,j)-dpcy*dyc(j);
                Pst=Psb(i);
                mdots=-Rho(i,j)*dx(i)*USiny;
                bbsy=-mdots*USiny+mbus(i)*(dx(i)/(dyc(j)))*USiny+Psb(i)*dx(i);
                ast=ast+mdots;
                elseif BSi==5
                elseif BSi==6
                elseif BSi==7
                mdots=-Rho(i,j)*dx(i)*Uy(i,j);
                Psb(i)=PSin;
                Pst=PSin;
                bbsy=PSout*dx(i);
                elseif BSi==8
                Psb(i)=P(i,j)-dpcy(i,j)*dyc(j);
                Pst=Psb(i);
                mdots=-Rho(i,j)*dx(i)*Uy(i,j);
                bbsy=-mdots*dux(i,j)*dyc(j)+Psb(i)*dx(i);
                end
                
            else
                dpcx(i,j)=(P(i+1,j)-P(i-1,j))/(dxc(i)+dxc(i+1));
                dpcy(i,j)=(P(i,j+1)-P(i,j-1))/(dyc(j)+dyc(j+1));
                dUys(i,j)=(Uy(i,j)-Uy(i,j-1))/dyc(j);
                as(i,j)=-mus(i,j)*(dx(i)/dyc(j));
                ast=as(i,j);
                Pnt=Pn(i,j);
                Pet=Pe(i,j);
                Pwt=Pw(i,j);
            end
            if j==m-1 % North
                    dpcy(i,j)=-(P(i,j-1)-PbN(i))/(dyc(j)+dyc(j+1));
                    dUyn(i,j)=-(Uy(i,j)-UbNy(i))/dyc(j+1);
                    if WNi==1
                    abx=-mbun(i)*(dx(i)/(dyc(j+1)));
                    Pnb(i)=P(i,j)+dpcy(i,j)*dyc(j+1);
                    Pnt=Pnb(i);
                    bbny=Pnb(i)*dx(i)-abx*UNwally;
                    ant=abx;
                    elseif WNi==2
                    Pnb(i)=P(i,j)+dpcy(i,j)*dyc(j+1);
                    bbny=Pnb(i)*dx(i);

                    elseif BNi==4
                    Pnb(i)=P(i,j)+dpcy*dyc(j+1);
                    Pnt=Pnb(i);
                    mdotn=Rho(i,j)*dx(i)*UNiny;
                    bbny=-mdotn*UNiny+mbun(j)*(dx(i)/(dyc(j+1)))*UNiny+Pnb(i)*dx(i);
                    ant=ant+mdotn;
                    elseif BNi==5
                    elseif BNi==6
                    elseif BNi==7
                    mdotn=Rho(i,j)*dx(i)*Uy(i,j);
                    Pnb(i)=PNin;
                    bbny=PNout*dx(i);
                    elseif BNi==8
                    Pnb(i)=P(i,j)+dpcy(i,j)*dyc(j+1);
                    Pnt=Pnb(i);
                    mdotn=-Rho(i,j)*dx(i)*Uy(i,j);
                    bbny=-mdotn*dux(i,j)*dyc(j+1)+Pnb(i)*dx(i);
                    end
                    
            else
                   dpcx(i,j)=(P(i+1,j)-P(i-1,j))/(dxc(i)+dxc(i+1));
                   dpcy(i,j)=(P(i,j+1)-P(i,j-1))/(dyc(j)+dyc(j+1));
                   dUyn(i,j)=(Uy(i,j+1)-Uy(i,j))/dyc(j+1); 
                   an(i,j)=-mun(i,j)*(dx(i)/dyc(j+1));
                   ant=an(i,j);
                   Pst=Ps(i,j);
                   Pet=Pe(i,j);
                   Pwt=Pw(i,j);
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
                    
                else
                    if Ux>0
                    if i==2
                        Uyuw=UWiny;
                    else
                        Uyuw=Uy(i-2,j);
                    end
                    Uydw=Uy(i,j);
                    Uycw=Uy(i-1,j);
                    awc=mdotw;
                    aw(i,j)=aw(i,j)+awc;
                    awt=aw(i,j);
                    else
                        if i==n-1
                            Uyuw=UEiny;
                        else
                            Uyuw=Uy(i+1,j);
                        end
                        Uydw=Uy(i-1,j);
                        Uycw=Uy(i,j);
                    end
                    Qw=Uydw-Uyuw;
                    if Qw~=0
                        Uydimw=(Uxcw-Uxuw)/Qw;
                    else
                        Uydimw=0;
                    end
                    Uyfdimw=Smart(Uydimw);
                    Uyfw=Uyfdimw*Qw+Uyuw;
                    bdcw=mdotw*(Uyfw-Uycw);
                end
                
                if j==1 %South
                    
                else
                    if Uy>0
                    if j==2
                        Uyus=USiny;
                    else
                        Uyus=Uy(i,j-2);
                    end
                    Uyds=Uy(i,j);
                    Uycs=Uy(i,j-1);
                    asc=mdots;
                    as(i,j)=as(i,j)+asc;
                    ast=as(i,j);
                    else
                        if j==m-1
                            Uyus=UNiny;
                        else
                            Uyus=Uy(i,j+1);
                        end
                        
                        Uyds=Uy(i,j-1);
                        Uycs=Uy(i,j);
                    end
                    Qs=Uyds-Uyus;
                    if Qs~=0
                        Uycdims=(Uycs-Uyus)/Qs;
                    else
                        Uycdims=0;
                    end
                    Uyfdims=Smart(Uycdims);
                    Uyfs=Uyfdims*Qs+Uyus;
                    bdcs=mdots*(Uyfs-Uycs);
                end
                
                if i==n-1 %East
                    
                else
                    if Ux<0
                        aec=mdote;
                        ae(i,j)=ae(i,j)+aec;
                        aet=ae(i,j);
                        if i==n-2
                            Uyue=UEiny;
                        else
                            Uyue=Uy(i+2,j);
                        end
                        Uyde=Uy(i,j);
                        Uyce=Uy(i+1,j);
                    else
                      
                        if i==1
                            Uyue=UWiny;
                        else
                            Uyue=Uy(i-1,j);
                        end
                        Uyde=Uy(i+1,j);
                        Uyce=Uy(i,j);
                    end
                    Qe=Uyde-Uyue;
                    if Qe~=0
                        Uycdime=(Uyce-Uyue)/Qe;
                    else
                        Uycdime=0;
                    end
                    Uyfdime=Smart(Uycdime);
                    Uyfe=Uyfdime*Qe+Uyue;
                    bdce=mdote*(Uyfe-Uyce);
                end
                
                if j==m-1 %North
                    
                else
                    if Uy<0
                        anc=mdotn;
                        an(i,j)=an(i,j)+anc;
                        ant=an(i,j);
                        if n==n-2
                            Uyun=UNiny;
                        else
                            Uyun=Uy(i,j+1);
                        end
                        Uydn=Uy(i,j);
                        Uycn=Uy(i,j+1);
                    else
                        if j==1
                            Uyun=USiny;
                        else
                            Uyun=Uy(i,j-1);
                        end
                        Uydn=Uy(i,j+1);
                        Uycn=Uy(i,j);
                    end
                    Qn=Uydn-Uyun;
                    if Qn~=0
                        Uycdimn=(Uycn-Uyun)/Qn;
                    else
                        Uycdimn=0;
                    end
                    Uyfdimn=Smart(Uycdimn);
                    Uyfn=Uyfdimn*Qn+Uyun;
                    bdcn=mdotn*(Uyfn-Uycn);
                end
                
            bbctot=bbsc+bbnc+bbec+bbwc; %convection boundary
            bbtot=bbn+bbs+bbe+bbw; %diffusion boundary
            bbvisc=bbwx+bbex+bbny+bbsy;
            bptot= Pet*dx(i)+Pwt*dx(i)+Pst*dy(j)+Pnt*dy(j);%pressure term  
            bttot=dUyw*dy(j)+dUye*dy(j)+dUyn*dx(i)+dUys*dx(i); %transpose diffusion term
                    
                
            if ConvectionScheme=="Upwind"
                bdctot=0;
            elseif ConvectionScheme=="Smart"
            bdctot=bdce+bdcw+bdcn+bdcs;%deffered correction
            end
            
            if TransientScheme=="BE" 
                a0=Rho(i,j)*dx(i)*dy(j)/dt;
            elseif TransientScheme=="CN"
                a0=Rho(i,j)*dx(i)*dy(j)/(dt/2);
            else
                a0=0;
            end

                
            a(i,j)=-(ant+ast+aet+awt)+(mdote+mdotw+mdots+mdotn)+a0;
            b(i,j)=bbvisc-bdctot+bbctot+bbtot+bptot+bttot+a0*Uyoldtrans(i,j)+Fb*dx(i)*dy(j);
        end
    end
end