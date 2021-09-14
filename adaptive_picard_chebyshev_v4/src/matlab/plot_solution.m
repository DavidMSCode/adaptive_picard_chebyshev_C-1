% Codeto plot Trajectory and Hamiltonian
% Robyn Woollands

clear
close all
clc

mu = 3.986004418e5;



cd ../
load output.txt
cd matlab
% load rk1210Out.txt
% output = rk1210Out;
Period = 5059;

figure(1)
plot(output(:,1)./Period,output(:,2:4),'o-','Linewidth',2)
title('Position')
xlabel('Time (s)')
ylabel('Position (km)')
legend('x','y','z')

figure(2)
plot(output(:,1)./Period,output(:,5:7),'o-','Linewidth',2)
title('Velocity')
xlabel('Time (s)')
ylabel('Velocity (km/s)')
legend('vx','vy','vz')

figure(3)
semilogy(output(:,1)./Period,output(:,8),'ro-')
grid on
title('Hamiltonian')
xlabel('Time (Orbits)')
ylabel('| E - E_o | / | E_o |')
set(gca, 'FontName', 'Helvetica','FontSize',16)

%%
l = length(output);
for i=1:l
    elts(i,:) = cspice_oscelt(output(i,2:7)',output(i,1),mu);
end
%%
ts = output(:,1);
ox = output(:,2);
oy = output(:,3);
fg = figure(4);
cla

pbaspect([1,1,1])
daspect([1 1 1])
drawnow
 %plot earth equator
theta = 0:pi/180:2*pi;
req = 6378.1;
xe = req.*cos(theta);
ye = req.*sin(theta);
plot(xe,ye,'g');
hold on

%plot individual orbits
q = elts(1,1);
e = elts(1,2);
a = q/(1-e);
P = 2*pi*sqrt(a^3/mu);
T1 = P;
i0 = 1;
i1 = find(ts<=T1,1,'last');
per = scatter(0,0,150,'r','filled');
apo = scatter(0,0,150,'b','filled');
count=1;

%Format Plot
legend(["Earth","Periapsis","Apoapsis"],'Autoupdate','off')
xlim([-7000 7000])
ylim([-7000 7000])
fg.Position = [100 100 1200 1200];
grid on
alts=[];
ts_per = [];
i=1;
while (ts(i1)-ts(i0))>=P-100
    %get data for current orbit
    xi = ox(i0:i1);         
    yi = oy(i0:i1);
    to = ts(i0:i1);
    r = sqrt(xi.^2+yi.^2);
    %find periapsis and apoapasis index
    i_per = find(r==min(r),1);
    i_apo = find(r==max(r),1);
    %plot orbit
    ax = plot(ox(i0:i1),oy(i0:i1),'k','linewidth',2);
    %move apoapsis and periapsis
    per.XData = xi(i_per);
    per.YData = yi(i_per);
    apo.XData = xi(i_apo);
    apo.YData = yi(i_apo);
    alts(i) = r(i_per)-req;
    ts_per(i) = to(i_per);
    uistack(per,'top')
    uistack(apo,'top')
    %change label
    title(sprintf('Frame # %i',count))
    %draw figure and wait for 60 fps
    axis equal
    drawnow
    pause(1/60)
%Setup next loop    
    i0 = i1;
    q = elts(i1,1);
    e = elts(i1,2);
    a = q/(1-e);
    P = 2*pi*sqrt(a^3/mu);       %estimate period for next orbit
    T1 = T1+P;
    i1 = find(ts<=T1,1,'last');                 %get index of end of orbit
    count=count+1;
    i=i+1;
%Make previous orbit transparent
    if i0<length(ts) && (ts(i1)-ts(i0))>=P-100
        ax.Color='b';
        ax.Color(4)=0.02;
        ax.LineWidth=1;
    end
end
as = elts(:,1)./(1-elts(:,2));
Ts = sqrt(as.^3./mu)*2*pi;

req = 6378.1;
hold on


%%
figure(5)
plot(alts)
grid on

 figure(6)
 plot(xe,ye,'g');
hold on
 plot(output(:,2),output(:,3))