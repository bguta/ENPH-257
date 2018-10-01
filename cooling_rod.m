close all;


L = 0.3;
dx = 0.01;
dt = 0.00001;

N = L / dx;

k = 250.0;  % W / m / K
c = 910.0 ; % J / kg / K
p = 2712.0 ; % kg / m^3

r = 0.01 ;% m radius 
Pin = 10  ;% W  power in
T_amb = 20 + 272.15 ; % ambient temp
k_c = 5  ;% W/m^2/K convection constant
epi = 1 ;
sigma = 5.67 * 10^(-8) ;  % W/m^2/K (stefan-Boltzmann constant)

C = k / (c * p);

x = linspace(0,L,N);
u = repelem([20.0 + 272.15],N-1);
T = [50.0 + 272.15, u];

t = 0;

while T(1) > 20.1
    if t >= 30000
        disp(t);
    end
for j = 1:1000
for i = 2:(N-1)
    T(i) = T(i) + C * dt * (T(i-1) - 2*T(i) +T(i+1) )/(dx ^ 2);
    
    T(i) = T(i) - 2 * dt * k_c * (T(i) - T_amb) / (c * p * r); %convection loss
    
    T(i) = T(i) - 2 * dt * epi * sigma * (T(i) ^(4) - T_amb^(4)) / (c * p * r) ; % radiation loss
    
    
end
T(1) = T(1) - C * dt * (T(1) - T(2))/ ( dx^2); % heat transfer
T(1) = T(1) + Pin * dt /(c * p * pi * r^(2) * dx) ; % power gain
% T(1) = T(1) - dt * k_c * (T(1) - T_amb) / (c * p); % convection loss
% T(1) = T(1) - dt * epi * sigma * (T(1)^(4) - T_amb^(4)) / (c * p * r); % radiation loss

T(1) = T(1) -dt*(2/(c*r*p) + 1/(c*p*dx))* (epi*sigma*(T(1)^(4) - T_amb^(4)) + k_c * (T(1) - T_amb));

T(N) = T(N) + C * dt * (T(N-1) - T(N))/ (dx ^2); % gain of heat
% T(N) = T(N) - dt * k_c * (T(N) - T_amb) / (c * p); % convection loss
% T(N) = T(N) - dt * epi * sigma * (T(N)^(4) - T_amb^(4)) / (c * p * r); % radiation loss
T(N) = T(N) -dt*(2/(c*r*p) + 1/(c*p*dx))*(epi*sigma*(T(N)^(4) - T_amb^(4)) + k_c * (T(N) - T_amb));
t = t + dt;
% if T(3) >= 25.0
%     disp("done")
%     disp(t)
%     return
% end
end


if t >= 0
    disp(T(1) - 272.15)
    q = plot(x,T,'b-');
    axis([0 0.3 10+272.15 60+272.15])
    title('Temperature of an Aluminum rod as a function of time and distance')
    xlabel('Distance (m)','FontSize',10);
    ylabel('Temperature (degrees C)','FontSize',10);
    pause(0.001)
end
end



