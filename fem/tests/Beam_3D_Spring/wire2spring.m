load spring/mesh.nodes
t=mesh(:,3);
s=(2+10*pi)*t-1;
x=cos(s);
y=sin(s);
z=0.1*s/pi;
dz=1.0-10*pi;

s1=0.0;
msk=find(s<s1);
x(msk)=cos(s1);
y(msk)=sin(s1);
z(msk)=s(msk);

s2=10*pi;
msk=find(s>s2);
x(msk)=cos(2*s2);
y(msk)=sin(2*s2);
z(msk)=s(msk)+dz;

mesh(:,3)=x;
mesh(:,4)=y;
mesh(:,5)=z;

n=size(s);


file_id = fopen('spring/mesh.nodes', 'w');
for i=1:n	
  fprintf(file_id,"%d %d %f %f %f\n",i,-1,x(i),y(i),z(i));
endfor
fclose(file_id);
