syms t;% parameters
L = 0.01975; % helmholtz coils 2 parallel - QUIVER DONE
N = 320;
r = 0.067;
I = 0.41;
B1 = [0,0,0];

dt = 200; %dt, the fineness of the coil readings
scale = 1;

%helix - finding dl
%add loop here to go through coil - include entire code

x1 = L*t/(2*pi*N);  %definition of the helix (solenoid)
y1 = r*cos(t);
z1 = r*sin(t);

dx1 = diff(x1); %derivatives wrt t
dy1 = diff(y1);
dz1 = diff(z1);

x2 = L*t/(2*pi*N)+0.067;
y2 = r*cos(t);
z2 = r*sin(t);

dx2 = diff(x2);
dy2 = diff(y2);
dz2 = diff(z2);

temp = -0.5:0.05:0.5;

%tail vector point
Bans = [size(temp),3];
%origin vector point
org = [size(temp),3];

ind = 0;

for dirx = -0.3:0.1:0.3
    for diry = -0.3:0.1:0.3
        for dirz = -0.3:0.1:0.3
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

hold on;

t = linspace(0,640*pi,100000);

x1 = L*t/(2*pi*N);  %definition of the helix (solenoid)
y1 = r*cos(t);
z1 = r*sin(t);

x2 = L*t/(2*pi*N)+0.06;
y2 = r*cos(t);
z2 = r*sin(t);

plot3(x1,y1,z1,x2,y2,z2)

axis equal;

quiver3(org(:,1),org(:,2),org(:,3),Bans(:,1),Bans(:,2),Bans(:,3))

grid on;

hold off;