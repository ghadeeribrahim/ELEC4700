function AddHCPAtomicArray(LAtoms,WAtoms,X0,Y0,VX0,VY0,RotAng,InitDist,Temp,Type)
global C
global x y AtomSpacing
global nAtoms % MinX MaxX MinY MaxY
global AtomType Vx Vy Mass0 Mass1

if Type == 0
    Mass = Mass0;
else
    Mass = Mass1;
end

L = ((LAtoms - 1) + 0.5) * AtomSpacing;
W = (WAtoms - 1) * sqrt(3) / 2 * AtomSpacing;

p =1;
% n = 4;
for i=1:LAtoms
    for j = 1:WAtoms
        y(end + 1) = (sqrt(3)*(j - 1) * AtomSpacing) / 2.0;
        if rem(j, 2) == 1
            x(end + 1) = (i - 1) * AtomSpacing;
        else
            x(end + 1) = AtomSpacing / 2.0 + (i - 1) * AtomSpacing;
        end
        p = p + 1;
    end
end

x = x - L/2;
y = y - W/2;

numAtoms = p-1;

x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + X0;
y(nAtoms + 1:nAtoms + numAtoms) = y(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + Y0;

AtomType(nAtoms + 1:nAtoms + numAtoms) = Type;

if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb*Temp/Mass);

    Vx(nAtoms + 1:nAtoms + numAtoms) = std0*randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0*randn(1, numAtoms);
end

Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;

end
function AddHCPAtomicBlob(LAtoms, X0, Y0, VX0, VY0, RotAng, InitDist, Temp, Type)
global C
global x y AtomSpacing
global nAtoms % MinX MaxX MinY MaxY
global AtomType Vx Vy Mass0 Mass1

if Type == 0
    Mass = Mass0;
else
    Mass = Mass1;
end

% L = ((LAtoms-1)+0.5)*AtomSpacing;
% W = (LAtoms-1)*sqrt(3)/2*AtomSpacing;

p = 1;

for j=1:LAtoms + 1
    for i = 1:2*LAtoms - 1 - j
        y(end + 1) = (sqrt(3) * ( j - 1) * AtomSpacing) / 2.0;
        x(end + 1) = (-(2 * LAtoms - j - 2) * AtomSpacing) / 2.0 + i * AtomSpacing;
        p = p + 1;
        if j ~= 1
            y(end + 1) = -y(end);
            x(end + 1) = x(end);
            p = p + 1;
        end
    end
end

x = x - (max(x) + min(x)) / 2;

numAtoms = p-1;

x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms)-0.5)*AtomSpacing*InitDist + X0*AtomSpacing;
y(nAtoms + 1:nAtoms + numAtoms) = y(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms)-0.5)*AtomSpacing*InitDist + Y0*AtomSpacing;

AtomType(nAtoms + 1:nAtoms + numAtoms) = Type;

if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb*Temp/Mass);

    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
end

Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;

end
function [ output_args ] = AddParticleStream(num, x0, y0, PartAng, Type, Ep, Seper)
global AtomSpacing x y AtomType Vx Vy Mass0 Mass1 nAtoms

if Type == 0
    Mass = Mass0;
else
    Mass = Mass1;
end

for p = 0:num - 1
    nAtoms = nAtoms + 1;
    x(nAtoms) = x0 * AtomSpacing - Seper * p * AtomSpacing * cos(PartAng);
    y(nAtoms) = y0 * AtomSpacing - Seper * p * AtomSpacing * sin(PartAng);
    AtomType(nAtoms) = Type;
end

V = sqrt(2 * Ep / Mass);

for p = 1:num
    Vx(nAtoms - num + p) = V * cos(PartAng);
    Vy(nAtoms - num + p) = V * sin(PartAng);
end

end
function GetForces(PhiCutoff,Epsilon,sigma)
global x y nAtoms Fx Fy Phi

% n = 1;
for i = 1:nAtoms
    Fx(i) = 0;
    Fy(i) = 0;
    Phi(i) = 0;
%     if i == 10
%         plot(x(i),y(i),'o','markers',48);
%         hold on
%     end

    for j = 1:nAtoms
        if i == j, continue; end
        dx = x(i) - x(j);
        dy = y(i) - y(j);
        r = sqrt(dx^2 + dy^2);

        if r > PhiCutoff, continue, end

        [aPhi dPhidr] = LJPot(r, Epsilon, sigma);

        ang = atan2(dy, dx);
%         dx =  r*cos(ang)
%         dy =  r*sin(ang)

        dFx = - dPhidr * cos(ang);
        dFy = - dPhidr * sin(ang);

%         if i == 10
%             plot(x(j),y(j),'ro','markers',48);
%         end

        Phi(i) = Phi(i) + aPhi;
        Fx(i) = Fx(i) + dFx;
        Fy(i) = Fy(i) + dFy;
    end
end

hold off

end


    