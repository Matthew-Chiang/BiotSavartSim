syms t;% parameters
syms thetha;
L = 0.01975; % helmholtz coils 2 parallel
N = 320;
r = 0.067;
I = 0.41;
B1 = [0,0,0];

dt = 500; %dt, the fineness of the coil readings
temp = -0.2:0.1:0.5 ;

scale = 0.5;

%helix - finding dl
%add loop here to go through coil - include entire code

x1 = L*t/(2*pi*N);  %definition of the helix (solenoid)
y1 = r*cos(t);
z1 = r*sin(t);

x2 = L*t/(2*pi*N);
y2 = r*cos(t);
z2 = r*sin(t);


Ry = [cos(thetha) 0 sin(thetha) ; 0 1 0 ; -sin(thetha) 0 cos(thetha)];
Ryin = inv(Ry);

Rz = [cos(thetha) -sin(thetha) 0 ; sin(thetha) cos(thetha) 0 ; 0 0 1];
Rzin = inv(Rz);

Rinit1 = [x1; y1 ;  z1];
Rinit2 = [x2 ;y2 ;z2];

Rfin1 = [];
Rfin2 = [];

Rtemp1 = Rz*Rinit1;
Rfin1 = Ry*Rtemp1;

Rfin2 = Rzin*Rinit2;

x1 = subs(Rfin1(1),thetha,pi/6);
y1 = subs(Rfin1(2),thetha,pi/6);
z1 = subs(Rfin1(3),thetha,pi/6);

x2 = subs(Rfin2(1),thetha,pi/6)+0.2;
y2 = subs(Rfin2(2),thetha,pi/6);
z2 = subs(Rfin2(3),thetha,pi/6);

dx1 = diff(x1); %derivatives wrt t
dy1 = diff(y1);
dz1 = diff(z1);

dx2 = diff(x2);
dy2 = diff(y2);
dz2 = diff(z2);
%tail vector point
Bans = [size(temp),3];
%origin vector point
org = [size(temp),3];

ind = 0;

for dirx = -0.2:0.1:0.5
    for diry = -0.4:0.1:0.3
        for dirz = -0.4:0.1:0.3
        ind = ind+1;
        ind
 % span of readings
        pos = [dirx diry dirz];
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
            rp1 = pos-l1; %rprime, distance from coil to point

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
            rp2 = pos-l2; %rprime, distance from coil to point

            C2 = cross(dl2,rp2);
	 % magnetic field at a certain point
            Btemp2 = 10^-7*I*C2/((norm(rp2))^3);
            B2 = B2+ Btemp2; %add all the magnetic fields due to each point together
            end%end loop here

        BT = B1+B2;
        
        org(ind,1) = dirx;
        org(ind,2) = diry;
        org(ind,3) = dirz;
        
        Bans(ind,1) = BT(1);
        Bans(ind,2) = BT(2);
        Bans(ind,3) = BT(3);
        
        end
    end
end
%scaling
   
for i4 = 1:length(Bans(:, 1))
        Bans(i4, :) = Bans(i4, :)/norm(Bans(i4,:))*scale;
end



t = linspace(0,640*pi,100000);

x3 = subs(Rfin1(1),thetha,pi/6)
y3 = subs(Rfin1(2),thetha,pi/6)
z3 = subs(Rfin1(3),thetha,pi/6)

x4 = subs(Rfin2(1),thetha,pi/6)+0.2
y4 = subs(Rfin2(2),thetha,pi/6)
z4 = subs(Rfin2(3),thetha,pi/6)   


x3 = (67*sin(t))/2000 - (3^(1/2)*((67*cos(t))/2000 - (79*3^(1/2)*t)/(5120000*pi)))/2;
y3 = (79*t)/(5120000*pi) + (67*3^(1/2)*cos(t))/2000;
z3 = (67*cos(t))/4000 + (67*3^(1/2)*sin(t))/2000 - (79*3^(1/2)*t)/(10240000*pi);

x4 = (67*cos(t))/2000 + (79*3^(1/2)*t)/(5120000*pi) + 1/5;
y4 = (67*3^(1/2)*cos(t))/2000 - (79*t)/(5120000*pi);
z4 = (67*sin(t))/1000;

hold on;

xlabel('X Position (m) =>', 'FontSize', 12)
ylabel('<= Y Position (m)', 'FontSize', 12)
zlabel('Z Position (m) =>', 'FontSize', 12)

plot3(x3,y3,z3,x4,y4,z4, 'Color', 'b')

axis equal;

quiver3(org(:,1),org(:,2),org(:,3),Bans(:,1),Bans(:,2),Bans(:,3),'Color',[1 0.549 0])

grid on;

hold off;