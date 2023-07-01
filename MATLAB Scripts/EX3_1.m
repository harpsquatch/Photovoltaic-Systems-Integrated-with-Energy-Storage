clc;clear;close all
%% Input data
f=logspace(5,9,1000); % 100 KHz~1 GHz with 1000 points
omega=2*pi*f;
Vs=1;
Zs=10;
L_load=50/(2*pi*30e6); % Z=100+50*i @ 30 MHz
ZL=100+1i*omega*L_load;
len=1;
l=0.5e-6;
c=200e-12;
%% Evaluate the line parameters
Zc=sqrt(l/c);
v=1/sqrt(l*c);
beta=omega/v;
gamma0=1i*beta;

phi11=cosh(gamma0.*len);
phi12=-sinh(gamma0.*len).*Zc;
phi21=-sinh(gamma0.*len)./Zc;
phi22=phi11;

%% solve the voltages and currents
I0=((phi11-ZL.*phi21).*Vs)./...
    ((ZL+Zs).*phi11-phi12-ZL.*phi21.*Zs);
V0=Vs-Zs.*I0;
VL=phi11.*V0+phi12.*I0;
Zin=V0./I0;

%% plot voltage
figure(1)
fig1=semilogx(f, db(abs(V0)),'r'...
    ,f, db(abs(VL)),'b');
set(gca, 'Fontsize', 18, 'Linewidth',1);
set(fig1, 'Linewidth', 4);
xlabel('Frequency, [Hz]');
ylabel('|V(0)|&|V(L)|,[dBV]');
title('Line voltage');
legend('V(0)','V(L)');
grid on;

%% plot input impedance
figure(2)
fig2=semilogx(f, db(abs(Zin)),'g');
set(gca, 'Fontsize', 18, 'Linewidth',1);
set(fig2, 'Linewidth', 4);
xlabel('Frequency, [Hz]');
ylabel('|Zin(0),[dB\Omega]');
title('Input impedance');
grid on
legend('Zin(0)');

figure(3)
fig2=semilogx(f, angle(Zin)*180/pi,'g');
set(gca, 'Fontsize', 18, 'Linewidth',1);
set(fig2, 'Linewidth', 4);
xlabel('Frequency, [Hz]');
ylabel('Phase(Zin),deg');
title('Input impedance');
grid on
legend('Zin(0)');

%% 30 MHz
[~,f0]=min(abs(f-30e6));
VL_30M=VL(f0)
MVL_30M=abs(VL(f0))



