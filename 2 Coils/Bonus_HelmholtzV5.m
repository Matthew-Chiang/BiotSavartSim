syms t;% parameters
L = 0.01975; % helmholtz coils 2 parallel
N = 320;
r = 0.067;
I = 0.41;
B1 = [0,0,0];

dt = 500; %dt, the fineness of the coil readings
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

%tail vector point
finalx = [];
finaly = [];
finalz = [];

%origin vector point
origx = [];
origy = [];
origz = [];

indx = 0;
indy = 0;
indz = 0;

for dirx = -1:0.5:0.5
indx = indx+1;
indy = 0;
    for diry = -1:0.5:0.5
    indy = indy+1;
    indz = 0;
        for dirz = -1:0.5:0.5
        indz = indz+1;
 % span of readings
        pos = [dirx diry dirz];
        B1 = [0,0,0];
        B2 = [0,0,0];
            for a1 = 0:dt:2*pi*N %going through the helix for one point (dt)
            a1   %coil 1

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
            a2   %coil 2

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
        
        origx(indx,indy,indz) = dirx;%for quiver
        origy(indx,indy,indz) = diry;
        origz(indx,indy,indz) = dirz;
        
        finalx(indx,indy,indz) = BT(2);
        finaly(indx,indy,indz) = BT(1);
        finalz(indx,indy,indz) = BT(3);

        end
    end
end
%scaling

for i = 1:length(origx(1, :))
        origx(:, i) = origx(:, i)/norm(origx(:, i))*scale;
end
for i2 = 1:length(origy(1, :))
        origy(:, i2) = origy(:, i2)/norm(origy(:, i2))*scale;
end
for i3 = 1:length(origz(1, :))
        origz(:, i3) = origz(:, i3)/norm(origz(:, i3))*scale;
end
    
for i4 = 1:length(finalx(1, :))
        finalx(:, i4) = finalx(:, i4)/norm(finalx(:, i4))*scale;
end
for i5 = 1:length(finaly(1, :))
        finaly(:, i5) = finaly(:, i5)/norm(finaly(:, i5))*scale;
end
for i6 = 1:length(finalz(1, :))
        finalz(:, i6) = finalz(:, i6)/norm(finalz(:, i6))*scale;
end
    


quiver3(origx,origy,origz,finalx,finaly,finalz)