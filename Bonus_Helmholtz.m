syms t;
L = 0.1;
N = 300;
dt = 0.001;
r = 0.06;
rv = [1,1,1];
I = 1;
B = [0,0,0];

%helix - finding dl
%add loop here to go through coil - include entire code

x = L*t;
y = r*cos(2*pi*N*t);
z = r*sin(2*pi*N*t);

dx = diff(x);
dy = diff(y);
dz = diff(z);

%tval = (pi*r^2)*(fine)/(dt) %how many segments of dt are there in the coil
%tval = pi*2*r*L/dt  %updated tval
f = (dx^2+dy^2+dz^2)^1/2;
helixL = int(f,0,L)/dt;

for a = 0:dt:helixL

dx = diff(x);
dy = diff(y);
dz = diff(z);

dl = a*[dx dy dz];

%r finding rprime

l = [x,y,z];
rp = rv-l;

C = cross(dl,rp);
Btemp = (10^-7)*I*C/(norm(rp))^3;

end%end loop here

B = B+ Btemp





