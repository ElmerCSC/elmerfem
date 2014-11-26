alp = 0.7;
mulx=0.0;
muly=0.1;

tmin=140;
tmax=200;

% get the range of movement
for t=1:tmax
  tn=sprintf('t%d_n.dat',t);
  n=load(tn);
  
  ta=sprintf('t%d_a.dat',t);
  a=load(ta);
  if t==1 
    xmin = min(a(:,1));
    xmax = max(a(:,1));
    ymin = min(a(:,2));
    ymax = max(a(:,2));
  else
    xmin = min(min(a(:,1)),xmin);
    xmax = max(max(a(:,1)),xmax);
    ymin = min(min(a(:,2)),ymin);
    ymax = max(max(a(:,2)),ymax);
  end  

  for m=1:n
    tb=sprintf('t%d_b%d.dat',t,m);
    b=load(tb);
    xmin = min(min(b(:,1)),xmin);
    xmax = max(max(b(:,1)),xmax);
    ymin = min(min(b(:,2)),ymin);
    ymax = max(max(b(:,2)),ymax);
  end 
end 

xmin
xmax
ymin
ymax

dx=xmax-xmin;
dy=ymax-ymin;
cx=0.1;
cy=0.05;

v=[xmin-cx*dx,xmax+cx*dx,ymin-cy*dy,ymax+cy*dy];

for t=tmin:tmax

tn=sprintf('t%d_n.dat',t);
n=load(tn);
  
ta=sprintf('t%d_a.dat',t);
a=load(ta);
  
cla
hold off
h=patch(a(:,1),a(:,2),'b');
hold on

for m=1:n
  tb=sprintf('t%d_b%d.dat',t,m);
  b=load(tb);

  te=sprintf('t%d_e%d.dat',t,m);
  e=load(te);

 if m==1 
   col='m';
 elseif m==2
   col='r';
 elseif m==3
   col='y';
 elseif m==4
   col='g';
 elseif m==5
   col='c';
 end

 patch(b(:,1),b(:,2),col);
 for j=2:size(e,1)-1
	 plot([e(1,1) e(j:j+1,1)' e(1,1)],[e(1,2) e(j:j+1,2)' e(1,2)],'k','LineWidth',1.0)
 end 
end 
  
alpha(alp);
axis(v);
axis off;


set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0 0 19.2 10.8]);

outfile = sprintf('tri%d.png',t);
print('-dpng','-r100',outfile);

%outfile = sprintf('tri%d.jpg',t);
%print('-djpeg','-r100',outfile);

%outfile = sprintf('tri%d.tif',t);
%print('-dtiff','-r100',outfile);


end % t-loop
