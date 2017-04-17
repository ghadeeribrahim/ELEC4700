%Assignment 2 part 1 b. 
%finite defference method %L/W=3/2
W=20;
L=30;

x=round(linspace(0,L,30));
y=round(linspace(0,W,20)); 
V=zeros(length(x), length(y));
Vx=zeros(length(x),1);
V0=5;

%BCs
V(1,:)=V0;
V(end,:)=V0;
V(:,1)=0;
V(:,end)=0; 


for n=0:400
    for ydx= 2:length(y)-1
        for idx=2:length(x)-1
            if ydx == 1
                V(idx,ydx) = V(idx,ydx+1);
            elseif ydx == length(y)
                V(idx,ydx) = V(idx,ydx-1);
            else           
                V(idx,ydx)= 1/4*(V(idx+1,ydx)+V(idx-1,ydx)+V(idx,ydx+1)+V(idx,ydx-1));
                Vx=V(:,5);

            end
        end
    end
    figure(1)
   surface(V.');pause(0.01)
end
figure(2)

plot(Vx)

[xx,yy]=meshgrid(x,y);
figure(3)
mesh(xx,yy,V.')
