itime=min(find(t>0.36));
t(itime)

plot(x(:,itime),b(:,itime))
hold on
plot(x(:,itime),p(:,itime)/20,'r')
plot(x(:,itime),n(:,itime)/10,'g')

axis([0 1 -0.2 1.2])
