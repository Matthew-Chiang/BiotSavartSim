syms t;% parameters
syms thetha; 

L = 0.01975; % helmholtz coils 2 parallel
N = 320;
r = 0.067;
I = 0.510247933;
B1 = [0,0,0];

dt = 250; %dt, the fineness of the coil readings

%helix - finding dl
%add loop here to go through coil - include entire code

x1 = L*t/(2*pi*N);  %definition of the helix (solenoid)
y1 = r*cos(t);
z1 = r*sin(t);

dx1 = diff(x1); %derivatives wrt t
dy1 = diff(y1);
dz1 = diff(z1);

x2 = L*t/(2*pi*N);
y2 = r*cos(t);
z2 = r*sin(t);

dx2 = diff(x2);
dy2 = diff(y2);
dz2 = diff(z2);

Ry = [1 0 0 ; 0 cos(thetha) -sin(thetha); 0 sin(thetha) cos(thetha)];

Rinit1 = [x1; y1 ;  z1];

Rfin1 = [];

Rfin1 = Ry*Rinit1;

x1 = subs(Rfin1(1),thetha,74*pi/180);
y1 = subs(Rfin1(2),thetha,74*pi/180);
z1 = subs(Rfin1(3),thetha,74*pi/180);

final = []; %final matrix for plotting
final1 = [];
final2 = [];
b1 = 0;

for temprv1 = -0.5:0.03:0.6%if change this also change plotx
 % span of readings
    b1 = b1+1
    rv1 = [temprv1 0 0];
    B1 = [0,0,0];
    B2 = [0,0,0];
for a1 = 0:dt:2*pi*N %going through the helix for one point (dt)

subx1 = subs(x1,a1); %sub a into parametric function(for part 2) 
suby1 = subs(y1,a1);
subz1 = subs(z1,a1);

subdx1 = subs(dx1,a1); %sub a into parametric function (for dl)
subdy1 = subs(dy1,a1);
subdz1 = subs(dz1,a1);

dl1 = dt.*[subdx1 subdy1 subdz1];

%r finding rprime

l1 = [subx1,suby1,subz1];
rp1 = rv1-l1; %rprime, distance from coil to point

C1 = cross(dl1,rp1);
	 % magnetic field at a certain point
Btemp1 = 10^-7*I*C1/((norm(rp1))^3);
B1 = B1+ Btemp1; %add all the magnetic fields due to each point together
end%end loop here

for a2 = 0:dt:2*pi*N %going through the helix for one point (dt)

subx2 = subs(x2,a2); %sub a into parametric function(for part 2) 
suby2 = subs(y2,a2);
subz2 = subs(z2,a2);

subdx2 = subs(dx2,a2); %sub a into parametric function (for dl)
subdy2 = subs(dy2,a2);
subdz2 = subs(dz2,a2);

dl2 = dt.*[subdx2 subdy2 subdz2];

%r finding rprime

l2 = [subx2,suby2,subz2];
rp2 = rv1-l2; %rprime, distance from coil to point

C2 = cross(dl2,rp2);
	 % magnetic field at a certain point
Btemp2 = 10^-7*I*C2/((norm(rp2))^3);
B2 = B2+ Btemp2; %add all the magnetic fields due to each point together
end%end loop here
B1
B2
BT = B1+B2

final1(1,b1) = norm (B1);
final2(1,b1) = norm (B2);
final(1,b1) = norm (BT);%insert magnitude of magnetic field into the final array
end

final
plotx = [-0.5:0.03:0.6];
plot(plotx,final)




