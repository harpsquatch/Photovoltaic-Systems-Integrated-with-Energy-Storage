

calculpri = calcul.*pri %Price is in cents/KWh, calcul is KWh
calculprisingle = reshape(calculpri,[],1);
calculprisinglesum = sum(calculprisingle)*30/100 % this is going to be in cents (30 for months)


