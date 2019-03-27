syms t;% parameters
L = 0.1;
N = 300;
r = 0.06;
I = 1;
B = [0,0,0];

dt = 300; %dt, the fineness of the coil readings

%helix - finding dl
%add loop here to go through coil - include entire code

x = L*t/(2*pi*N);  %definition of the helix (solenoid)
y = r*cos(t);
z = r*sin(t);

dx = diff(x); %derivatives wrt t
dy = diff(y);
dz = diff(z);

%tval = (pi*r^2)*(fine)/(dt) %how many segments of dt are there in the coil
%tval = pi*2*r*L/dt  %updated tval
f = (dx^2+dy^2+dz^2)^1/2; %arc length of helix
helixL = int(f,0,2*pi*L); %divide into dt parts

final = []; %final matrix for plotting
b = 0;

for temprv = -0.2:0.5:0.2%if change this also change plotx
 % span of readings
    b = b+1
    rv = [temprv 0 0]
    B = [0,0,0];
for a = 0:dt:2*pi*N %going through the helix for one point (dt)

subx = subs(x,a) %sub a into parametric function(for part 2) 
suby = subs(y,a)
subz = subs(z,a)

subdx = subs(dx,a) %sub a into parametric function (for dl)
subdy = subs(dy,a)
subdz = subs(dz,a)

dl = dt.*[subdx subdy subdz]

%r finding rprime

l = [subx,suby,subz]
rp = rv-l %rprime, distance from coil to point

C = cross(dl,rp)
	 % magnetic field at a certain point
Btemp = 10^-7*I*C/((norm(rp))^3);
B = B+ Btemp; %add all the magnetic fields due to each point together
end%end loop here

final(1,b) = norm(B)%insert magnitude of magnetic field into the final array
end

final
plotx = [-0.2:0.5:0.2];
plot(plotx,final)




