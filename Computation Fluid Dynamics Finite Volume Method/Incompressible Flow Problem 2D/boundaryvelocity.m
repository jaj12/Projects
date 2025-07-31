function [UbWx,UbWy,UbEx,UbEy,UbSx,UbSy,UbNx,UbNy]= boundaryvelocity(UWwallx,UWwally,UEwallx,UEwally,UNwallx,UNwally,USwallx,USwally,BE,BW,BN,BS,WN,WE,WS,WW,Ux,Uy,UNiny,USiny,UNinx,USinx,UEinx,UWinx,n2i,s2i,e2i,w2i,n,m)

for i=1:n-1
    for j=1:m-1
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
           
                
            if WWi==9
                UbWx(j)=UWwallx;
                UbWy(j)=UWwally;
            elseif WWi==10
                UbWx(j)=Ux(i,j);
                UbWy(j)=Uy(i,j);
            elseif BWi==4
                UbWx(j)=UWinx;
                UbWy(j)=UWiny;
            elseif BWi==5
                
            elseif BWi==6
                
            elseif BWi==7
                UbWx(j)=Ux(i,j);
                UbWy(j)=Uy(i,j);
            elseif BWi==8
                UbWx(j)=Ux(i,j);
                UbWy(j)=Uy(i,j);
            end
        end
        
        if i==n-1
            
                
            if WEi==9%East
                UbEx(j)=UEwallx;
                UbEy(j)=UEwally;
            elseif WEi==10
                UbEx(j)=Ux(i,j);
                UbEy(j)=Uy(i,j);
            elseif BEi==4
                UbEx(j)=UEinx;
                UbEy(j)=UEiny;
            elseif BEi==5
                
            elseif BEi==6
                
            elseif BEi==7
                UbEx(j)=Ux(i,j);
                UbEy(j)=Uy(i,j);
            elseif BEi==8
                UbEx(j)=Ux(i,j);
                UbEy(j)=Uy(i,j);
            end
        end
        
        if j==1%South
            if WSi==9
                UbSx(i)=USwallx;
                UbSy(i)=USwally;
            elseif WSi==10
                UbSx(i)=Ux(i,j);
                UbSy(i)=Uy(i,j);
            elseif BSi==4
                UbSx(i)=USinx;
                UbSy(i)=USiny;
            elseif BSi==5
                
            elseif BSi==6
                
            elseif BSi==7
                UbSx(i)=Ux(i,j);
                UbSy(i)=Uy(i,j);
            elseif BSi==8
                UbSx(i)=Ux(i,j);
                UbSy(i)=Uy(i,j);
            end
            
        end
        
        if j==m-1 %North
            if WNi==9
                UbNx(i)=UNwallx;
                UbNy(i)=UNwally;
            elseif WNi==10
                UbNx(i)=Ux(i,j);
                UbNy(i)=Uy(i,j);
            elseif BNi==4
                UbNx(i)=UNinx;
                UbNy(i)=UNiny;
            elseif BNi==5
                
            elseif BNi==6
                
            elseif BNi==7
                UbNx(i)=Ux(i,j);
                UbNy(i)=Uy(i,j);
            elseif BNi==8
                UbNx(i)=Ux(i,j);
                UbNy(i)=Uy(i,j);
            end
        end
    end
end