% Data provided:
Vs = 1; %source Voltage
f = logspace(8,9,1000); %frequency interval
w = 2*pi*f; %omega
l = 0.3372797e-6; % per unit inductance
c = 62.649e-12; %per unit Capacitance
L = 5; %length of transmission line
%Line parameters:
Zc = sqrt(l/c); % Chracteristic impedence
v = 1/sqrt(l*c); %Velocity
beta = w./v; %phase
gamma= j*beta; %propagation constant
%Chain Parameters
phi11 = cosh(gamma.*L);
phi21= -sinh(gamma.*L)./Zc;
phi12= Zc.*-sinh(gamma.*L);
phi22= phi11;
%Ploting line voltages:
for i=1:1:3
figure
if i==1 % Case A
ZS = 50;
ZL = 50;
end
if i==2 % Case B
ZL = 75;
ZS = 75;
end
if i ==3 % Case C
ZL = 1000;
ZS = 1000;
end
% Calculations
I0 = ((phi11 - ZL.* phi21) .* Vs)./((ZL+ZS).* phi11 - phi12 - ZL .* phi21 .* ZS); V0 = Vs - ZS .* I0;
VL = phi11 .* V0 + phi12 .* I0;
RC=(ZL-Zc)/(ZL+Zc);%reflection coefficient
VSWR=(1+RC)/(1-RC); % Voltage standing wave ratio
fig=semilogx( f, 20*log10(abs(V0)) ,'k', f, 20*log10(abs(VL)) ,'r');
set(gca, 'FontSize',12,'LineWidth',1);
set(fig,'LineWidth', 2);
xlabel('Frequency, [Hz]');
ylabel('|V(0)| & |V(L)|, [dBV]');
title(['Case: ' num2str(i) ' Line voltages VSWR ' num2str(VSWR)]);
end