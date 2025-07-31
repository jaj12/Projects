function [TbN,TbS,TbE,TbW]=borderTemp(Tc,TbN,TbS,TbE,TbW,BW,BE,BN,BS,QbW,QbE,QbS,QbN,hbW,hbE,hbS,hbN,TinfW,TinfE,TinfS,TinfN,dx,dy,dxc,dyc,TWin,TEin,TNin,TSin)
    kbn= cond(TbN);
    kbs= cond(TbS);
    kbe= cond(TbE);
    kbw= cond(TbW);
    if BW==2
        TbW=QbW.*dy./kbw+Tc(1,:);
    elseif BW==3 
        TbW=(hbW.*dy.*TinfW+kbw.*(dy/dxc(1)).*Tc(1,:))./(hbW.*dy+kbw.*(dy/dxc(1)));
    elseif BW==4 
        TbW=TWin*ones(length(TbW),1);
    elseif BW==5
        TbW=Tc(1,:);
    end
    
    if BE==2
        TbE=QbE.*dy./kbe+Tc(end,:);
    elseif BE==3
        TbE=(hbE.*dy.*TinfE+kbe.*(dy/dxc(end)).*Tc(end,:))./(hbE.*dy+kbe.*(dy/dxc(end)));
    elseif BE==4
        TbE=TEin*ones(length(TbE),1);
    elseif BE==5
        TbE=Tc(end,:);
    end
    
    if BN==2
        TbN=QbN.*dx./kbn+Tc(:,end);
    elseif BN==3
        TbN=(hbN.*dx.*TinfN+kbn.*(dx/dyc(end)).*Tc(:,end))./(hbN.*dx+kbn.*(dx/dyc(end)));
    elseif BN==4
        TbN=TNin*ones(length(TbN),1);
    elseif BN==5
        TbN=Tc(:,end);
    end
    
    if BS==2
        TbS=QbS.*dx./kbs+Tc(:,1);
    elseif BS==3
        TbS=(hbS.*dx.*TinfS+kbs.*(dx/dyc(1)).*Tc(:,1))./(hbS.*dx+kbs.*(dx/dyc(1)));
    elseif BS==4
        TbS=TSin*ones(length(TbS),1);
    elseif BS==5
        TbS=Tc(:,1);
    end
    
end