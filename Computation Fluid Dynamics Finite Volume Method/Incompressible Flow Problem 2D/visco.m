function [y]= visco(x)
%conduction relationship to temperature given
%y=50.+x./10;
[n,m]=size(x);
y=5e-5*ones(n,m);
end