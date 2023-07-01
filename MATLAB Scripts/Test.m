calculprilowda = PSSLowda.*Prilowda %Price is in cents/KWh, calcul is KWh
calculprisingle = reshape(calculprilowda,[],1);
calculprisinglesumlow = sum(calculprisingle)*30/100 % this is going to be in cents so converting, *30 is for months

calculprilowda = PSSLowda.*PriESS %Price is in cents/KWh, calcul is KWh
calculprisingle = reshape(calculprilowda,[],1);
calculprisinglesumlowESS = sum(calculprisingle)*30/100 % this is going to be in cents so converting, *30 is for months


