
%Only PV 

Sexy = (P_pv(irr,k(tech_sel,:),y(tech_sel,find(self_cons_thres == self_cons_comp(1))),rend_PV)./1e3)
calcu = Sexy - p_load/1e3 

calcul = reshape(calcu/4,[57,12]) %this is the 57X12 ki matrix with difference in values  and is in KWh 

pri = zeros(size(calcul)) %this is in cent/KWh 



for ii = 1:57
    for jj = 1:12
        if calcul(ii,jj)>0
            pri(ii,jj) = 12.5             %this loop creates the PV window prices in cents/KWh
        else 
            pri(ii,jj) = 30
        end
    end
end 


loadprof = residential*(P_load/1e5) %Here im taking the load profile and multiplying it with P_load/(10^3 (Watt- KW) and 10^2(for percent) in KW

for ii = 21:76
        loadprof(ii) = 0  %This loop removes all the values from 5 to 7 
end 
 

PSSLowda = zeros(96,12) %Creating 96 X12 matrix for power, idea is to create a zero matrix for each representative day in each month 
Prilowda = zeros(96,12) % Doing the same for price

for ii = 1:size(PSSLowda,1)
    for jj = 1:size(PSSLowda,2)
        PSSLowda(ii,jj) = -loadprof(ii)/4 %with this loop im basically pasting everything in -ve (0(pv) - load) in this year wali matrix in KWh
    end 
end 

for aa = 1:size(calcul,1)
    for bb = 1:size(calcul,2) 
        PSSLowda(20+aa,bb) = calcul(aa,bb) %with this loop im pasting the window difference in KWh
    end 
end 

loadprofyr = repmat(residential*(P_load/1e5),1,12) 

for aa = 1:size(PSSLowda,1)
    for bb = 1:size(PSSLowda,2) 
        if PSSLowda(aa,bb) == 0
            PSSLowda(aa,bb) = -loadprofyr(aa,bb)/4   %with this loop im filling out the missing value in KWh
        end 
    end 
end

for aa = 1:size(calcul,1)
    for bb = 1:size(calcul,2) 
        Prilowda(20+aa,bb) = pri(aa,bb)    %with this loop im filling the PV window prices in big matrix in cents/KWh
    end 
end 

for aa = 1:size(Prilowda,1)
    for bb = 1:size(Prilowda,2) 
        if Prilowda(aa,bb) == 0
            Prilowda(aa,bb) = 30    %filling missing values in the price matrix 
        end 
    end 
end

calculprilowda = PSSLowda.*Prilowda %Price is in cents/KWh, calcul is KWh
calculprisingle = reshape(calculprilowda,[],1);
calculprisinglesumlow = sum(calculprisingle)*30/100 % this is going to be in cents so converting, *30 is for months
