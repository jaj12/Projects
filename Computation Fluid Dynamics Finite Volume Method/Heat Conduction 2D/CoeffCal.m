function [a,an,as,aw,ae,b]=CoeffCal(dxc,dyc,dy,dx,Qg,BE,BW,BN,BS,hbE,hbW,hbN,hbS,TinfE,TinfW,TinfN,TinfS,TbE,TbW,TbN,TbS,QbE,QbW,QbN,QbS,kn,ks,ke,kw)
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

                else
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

                else
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
                else

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

                    else
                    RbN=(hbN*kbn(i)/dyc(j+1))/(hbN+kbn(i)/dyc(j+1));
                    ab=-RbN*dx(i);
                    bbn=-ab*TinfN;
                    ant=ab;
                    end
                    
            else
                   an(i,j)=-kn(i,j)*(dx(i)/dyc(j+1));
                   ant=an(i,j);
            end
            
            a(i,j)=-(ant+ast+aet+awt);
            b(i,j)=bbn+bbs+bbe+bbw+Qg*dx(i)*dy(j);
        end
    end
end