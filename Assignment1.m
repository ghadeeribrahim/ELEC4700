%Constants 

n=1000; %number of electrons
m=0.26*9.1*10e-31; %electron mass
K=1.38*10e-23;  %Boltzmann constant 
T=300; %Temperature

%%Setting up the boundaries 
L=200*1e-9; 
W=100*1e-9;

%%Calculation of thermal velocity 
vth=sqrt(2*K*T/m); 
%%
% 
%  PREFORMATTED
%  TEXT
% 
%%vth= 1.8707e+05

%%
%assigning particles random position and velocity in the x-y plane 
xpos= rand(1,n)*L;
ypos= rand(1,n)*W;
xvel= randn(1,n)*vth;
yvel= randn(1,n)*vth;

%plot of the initial positions. 
figure(1)
plot(xpos,ypos,'o')
hold on

%Ensuring that we have a normal distribution of velocities by plotting a
%histogram 
histogram(sqrt(xvel.*xvel+yvel.*yvel),50)
hold on 

%setting up the for loop.     
dt = L/100/vth;
nsteps = 1000; 
Tcalc= zeros (n,1);

% for i = 1:nsteps
%     
% %    add new position to old position and assigning it to be the final
% %    position
%     x = xpos + dt*xvel;
%     y = ypos + dt*yvel;
%     
%     xshift= x>200e-9;
%     x(xshift)=x(xshift)-200e-9;
%     ux= x<0 ;
%     x(ux)= x(ux)+200e-9;  
%  
%     yshift=y>100e-9;
%     yvel(yshift)=-yvel(yshift); 
%     uy= y<0; 
%     yvel(uy)= -yvel(uy);
% 
%     plot(x,y,'o');
%     title ('Plot of electrons without scatter') 
%     axis([0 L 0 W])
%     pause(0.1)
%     
%    xpos = x;   
%    ypos = y;
% end
% hold off 
% 
% for i=1:nsteps
%     
%     x = xpos + dt*xvel;
%     y = ypos + dt*yvel;
%     
%     xshift= x>200e-9;
%     x(xshift)=x(xshift)-200e-9;
%     ux= x<0 ;
%     x(ux)= x(ux)+200e-9;  
%  
%     yshift=y>100e-9;
%     yvel(yshift)=-yvel(yshift); 
%     uy= y<0; 
%     yvel(uy)= -yvel(uy);
%       plot(x,y,'o'); 
%     axis([0 L 0 W])
%     pause(0.1)
%     hold on
%   
%   pscat= (1-exp(-dt/0.2e-12)); 
%     scat=rand(n,1); 
%     iscat=scat<pscat;
%     nones=nnz(iscat);
%     xvel(iscat)=sqrt(2.*K.*T./m).*randn(nones,1).*2;
%     yvel(iscat)=sqrt(2.*K.*T/m).*randn(nones,1).*2;
% 
%     figure (2)
%     plot(x,y,'o'); 
%     title ('Plot of Electrons with Scatter')
%     axis([0 L 0 W])
%     pause(0.1)
    
    vthsq=(xvel.^2)+(yvel.^2);
    vthsqm=vthsq.*m;
    Tcalc= (0.5.*vthsqm)./K;
    
    figure(3)
    plot(Tcalc)
    title('temperature plot') 
  
%    xpos = x;   
%    ypos = y;
% 
% %    figure(4) 
% %    histogram(sqrt(xvel.*xvel+yvel.*yvel),50)
% %    title('velocity histogram') 
% %    
% end

%   
