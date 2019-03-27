t = linspace(0,640*pi,100000);

L = 0.01975;
N = 320;
r = 0.067;
I = 0.41;

x1 = L*t/(2*pi*N);  %definition of the helix (solenoid)
y1 = r*cos(t);
z1 = r*sin(t);

x2 = L*t/(2*pi*N)+0.06;
y2 = r*cos(t);
z2 = r*sin(t);

temprv1 = -0.3:0.005:0.3
rv1 = [temprv1 0 0]

plot3(x1,y1,z1,x2,y2,z2)
