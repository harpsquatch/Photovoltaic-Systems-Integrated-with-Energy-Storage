SOCLowda =  zeros(96,12)



%charging 
for jj = 1:size(PSSLowda,2)
    for ii = 1:size(PSSLowda,1)
        if PSSLowda(ii,jj) > 0
            SOCLowda(ii,jj) = SOCLowda(ii-1,jj)+ (0.25)*((PSSLowda(ii,jj)*4)/9.2939)
        end 
    end
end    

%charging 
for jj = 1:size(PSSLowda,2)
    for ii = 1:size(PSSLowda,1)     % in this case we are actually 
        if SOCLowda(ii,jj) > 1
            SOCLowda(ii,jj) = 1
        end 
    end
end 

%discharging
for jj = 1:12
    for ii = 1:95
        if SOCLowda(ii,jj) > 0 && SOCLowda(ii+1,jj) == 0
            SOCLowda(ii+1,jj) = SOCLowda(ii,jj)+ (0.25)*((PSSLowda(ii+1,jj)*4)/9.2939)
        end
    end
end           
                
SOCLowda(SOCLowda<0) = zeros(size(find(SOCLowda<0)));

PriESS = Prilowda

for aa = 1:size(SOCLowda,1)
    for bb = 1:size(SOCLowda,2) 
        if PSSLowda(aa,bb) > 0
            PriESS(aa,bb) = 0   %filling missing values in the price matrix 
        end 
    end 
end

calculprilowda = PSSLowda.*PriESS %Price is in cents/KWh, calcul is KWh
calculprisingle = reshape(calculprilowda,[],1);
calculprisinglesumlowESS = sum(calculprisingle)*30/100 % this is going to be in cents so converting, *30 is for months

