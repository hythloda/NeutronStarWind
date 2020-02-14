clear
load out
nx     =out(1,2);
nxp    =nx+1;
ll     =length(out(:,1));
time_no=ll/nxp
t      =out(1:nxp:ll,1);
for i=1:time_no
    ii     =(2:nxp)+(i-1)*nxp;
    x(:, i)=out(ii,1);
    b(:, i)=out(ii,2);
    p(:, i)=out(ii,3);
    v(:, i)=out(ii,4);
    e(:, i)=out(ii,5);
end
g=(1-b.^2).^(-1/2);
n=v.^(-1);

%dv     =diff(x.^3)*4*pi/3;
%dv     =[dv;zeros(1,time_no)];
%etot   =g.^2 .*(e+b.^2 .*p) .* dv*1.5d-3;
%mtot   =g.* n .* dv;
%sumetot=sum(etot,1);
%summtot=sum(mtot,1);
