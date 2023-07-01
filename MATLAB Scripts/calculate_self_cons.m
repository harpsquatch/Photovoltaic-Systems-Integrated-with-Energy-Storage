function [ y ] = calculate_self_cons ( x , w , z , k, eta )
% CALCULATE_SELF_CONS evaluates the self-consumption of power system with
% PV and load
% y = calculate_self_cons ( x , w , z , k, eta )
% 
% OUTPUT
%   y - is the self consumption
% INPUT
%   x - is the nominal power of the PV plant
%   w - is the load profile
%   z - is the irradiance
%   k - is the 3-by-1 vector of parameters for generating power from
%       irradiance
% eta - is the efficiency of the PV plant

P_pv = max(0,eta*(k(3).*z.^2 + k(2).*z + k(1)).*x);
%this calculates the the power by your PV with the given irradiation 

a = sum(P_pv);
%Calclautes your total PV, sums up all your PV data for each column, each
%column represents your month i guess, the sum is like total energy
%produced by PV in a day for that month. 


b = P_pv;


b((b - w) > 0) = w ( (b - w) > 0 );
%this checks the ppv - pload wala formula 
%if Pv output is greater than the load then Psc is load value 
% if pv output is less than the load then Psc is Pv value 

y = sum(b)/a;
%whatever it is, you get the value of self consumption and this should be
%0.8, which means pload/pv should be 0.8 or 80%, 80% to load, 20% to
%battery 
end