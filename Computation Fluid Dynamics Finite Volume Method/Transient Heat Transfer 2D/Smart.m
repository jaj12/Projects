function [Tdimf]=Smart(Tdimc)
        if Tdimc>=0 && Tdimc<=1/6
                Tdimf=3*Tdimc;
            elseif Tdimc>=1/6 && Tdimc<=7/10
                Tdimf=3/4*Tdimc+3/8;
            elseif Tdimc>=7/10 && Tdimc<=1
                Tdimf=1/3*Tdimc+2/3;
            else
                Tdimf=Tdimc;
        end

%         if Tdimc>=0 && Tdimc<=5/6
%                 Tdimf=3/4*Tdimc+3/8;
%             elseif Tdimc>=5/6 && Tdimc<=1
%                 Tdimf=1;
%             else
%                 Tdimf=Tdimc;
%         end

end