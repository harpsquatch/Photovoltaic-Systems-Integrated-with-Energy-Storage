%% Topic 2
% These are some functions, operations, and analyses that could be used
% when tackling Topic 2


% General constants
g = 9.806;
rho = 1;
eta = .9*.99*.99*.99;
% Run-of-the-river hydro power plant
Hm = 40;
q1 = [20 20 25 30 55 55 50 45 45 40 30 20];
% Reservoir hydro power plant
H = 650;
V_2_tot = 250e6;

%% powers and energies in MW and MWh
% $P = g\rho \eta q H $
P = @(x,y) g.*rho.*eta.*y.*x.*1e-3;
% $W = \frac{g\rho \eta V H}{3600} $
W = @(x,y) P(x,y)./3600;
% $q = \frac{P}{g \rho \eta H}$
P_inv = @(x,y) (x.*1e3)./(g.*rho.*eta.*y);
compute_V = @(x,y) P_inv(x,y).*3600;

% available energy in P2 [GWh]
W(V_2_tot,H).*1e-3;
% equivalent constant power for the whole year [MW]
P2_ave = W(V_2_tot,H)./8760;
% equivalent flow [m^3/s]
P_inv(W(V_2_tot,H)./8760,H);


% Load
wdl_pu = horzcat(.5*ones(1,7), .7, .85, ones(1,4), .85, .85, ones(1,4), .5*ones(1,5))';
phl_pu = horzcat(.5*ones(1,8), .75*ones(1,5), .5*ones(1,11))';
hl_pu = .3*ones(1,24)';
pu_prof = horzcat(wdl_pu,phl_pu,hl_pu);
Pn = 160;
% Time-of-use structure
year = 2020;
sd = ['01/01/2020';'06/01/2020';'05/04/2020';'25/04/2020';'01/05/2020';'02/06/2020';'15/08/2020';'01/11/2020';'08/12/2020';'25/12/2020';'26/12/2020'];
[fasce_h,fasce_q] = time_class(year,sd);
% a general workig day structure
wd_tou = fasce_h(1*24+1:2*24); % IDX_TOU = 1;
ph_tou = fasce_h(3*24+1:4*24); % IDX_TOU = 2;
h_tou = fasce_h(1:24); % IDX_TOU = 3;
idx_tou = horzcat(wd_tou,ph_tou,h_tou);
%%
% building structures to store data
Pg = zeros(3,24,12);
P1 = Pg;
P2 = Pg;
DeltaV = zeros(3,25,12);
ctrl_infl = Pg;
nat_infl = Pg;
no_days = zeros(3,1,12);
used_volume = zeros(3,12);
pondage_V = used_volume;
%%
start_day = datenum(strcat('1/1/',num2str(year)),'dd/mm/yyyy');
end_day = start_day + yeardays(year)-1;
days_in_year(:,1) = start_day:1:end_day;
dummy = weekday(days_in_year);
days_classification = ones(size(dummy));
days_classification(dummy == 1) = 3;
days_classification(dummy == 7) = 2;
idx_holidays = find(ismember(days_in_year,datenum(sd,'dd/mm/yyyy'))==1);
days_classification(find(idx_holidays == 1)) = 3;
idx_holidays(find(idx_holidays == 1)) = [];
days_classification(idx_holidays) = 3;
days_classification(idx_holidays-1) = 2;

avail_vol_2 = V_2_tot;
done_q = [];
FIGURES = false;
%%
for kk = 1:1:12,
    if ismember(q1(kk),done_q)
        idx = find(q1 == q1(kk));
        fprintf(1,"%s has the same inflow of %s \n",datestr(datenum(strcat('1/',num2str(kk),'/',num2str(year)),'dd/mm/yyyy'),'mmm'),...
            datestr(datenum(strcat('1/',num2str(idx(1)),'/',num2str(year)),'dd/mm/yyyy'),'mmm'));
        fprintf(1,"updating arrays and subtracting volume to P2, according to the number of days, and then skipping all calculations\n\n");
        for hh =1:1:3,
            P1(hh,:,kk) = P1(hh,:,idx(1));
            P2(hh,:,kk) = P2(hh,:,idx(1));
            Pg(hh,:,kk) = Pg(hh,:,idx(1));
            ctrl_infl(hh,:,kk) = ctrl_infl(hh,:,idx(1));
            nat_infl(hh,:,kk) = nat_infl(hh,:,idx(1));
            DeltaV(hh,:,kk) = DeltaV(hh,:,idx(1));
            pondage_V(hh,kk) = pondage_V(hh,idx(1));
            used_volume(hh,kk) = used_volume(hh,idx(1));
            avail_vol_2 = avail_vol_2 - used_volume(hh,kk).*no_days(hh,kk);
        end
        clear idx
    else
        done_q = [done_q q1(kk)];
        if FIGURES
            figure(4*(kk-1)+1);
            figure(4*(kk-1)+2);
            figure(4*(kk-1)+3);
            figure(4*(kk-1)+4);
        end
        for hh = 1:1:3
            no_days(hh,kk) = sum(days_classification(find(month(days_in_year)==kk))==hh);
            P_load = pu_prof(:,hh)'.*Pn;
            P_bl = min(P_load);
            idx_ps(1,:) = find(((idx_tou(:,hh) == hh) & (pu_prof(:,hh) == max(pu_prof(:,hh)))) == 1);
            peak_start = idx_ps(1);
            peak_end = idx_ps(end);
            peak_hours = [peak_start:1:peak_end];
            E_ph_tot = sum(pu_prof(peak_hours,hh)'.*Pn - P_bl.*ones(size(peak_hours)));
            E_ph_res = sum(pu_prof(idx_ps,hh)'.*Pn - P_bl.*ones(size(idx_ps)));
            dummy = {peak_hours,idx_ps};
            % run-of-the-river hydro
            % P1 = zeros(1,24);
            % total available energy in a day
            E1_av = P(q1(kk),Hm).*24;
            % if the energy is enough to cover the peak load, then do it, otherwise,
            % focus on Fhh & max load
            idx_p1 = dummy{(E1_av>=E_ph_res)*1 + ~(E1_av>=E_ph_res)*2};
            if E1_av > E_ph_tot
                P1(hh,idx_p1,kk) = E_ph_tot./length(idx_p1);
                P1(hh,~ismember([1:24],idx_p1),kk) = (E1_av - E_ph_tot)./24;
                P1(hh,idx_p1,kk) = P1(hh,idx_p1,kk) + (E1_av - E_ph_tot)./24;
            else
                P1(hh,idx_p1,kk) = E1_av./length(idx_p1);
            end
            % compute the profiles of the volumes
            % volumes of the controlled inflow
            ctrl_infl(hh,:,kk) = cumsum(compute_V(P1(hh,:,kk),Hm));
            % volumes of the natural inflow
            nat_infl(hh,:,kk) = cumsum(q1(kk).*ones(1,24).*3600);
            if FIGURES
                figure(4*(kk-1)+1);
                subplot(3,1,hh)
                hold on
                p(1) = plot([0:24],horzcat(0,nat_infl(hh,:,kk)).*1e-6,'LineWidth',2);
                p(2) = plot([0:24],horzcat(0,ctrl_infl(hh,:,kk)).*1e-6,'LineWidth',2);
                set(gca,'XLim',[0 24],'XTick',[0:4:24],'FontSize',11,'FontWeight','demi')
                xlabel('time [h]')
                ylabel('volume [Mm^3]')
                legend('natural inflow','controlled inflow','Location','best')
                grid on
                box on
                hold off
            end
            % difference between natural inflow straight line and controlled inflow
            % curve
            DeltaV(hh,:,kk) = horzcat(0,nat_infl(hh,:,kk)) - horzcat(0,ctrl_infl(hh,:,kk));
            % volume of the pondage in Mm^3
            if isempty((max(DeltaV(hh,DeltaV(hh,:,kk) > 0,kk)) + abs(min(DeltaV(hh,DeltaV(hh,:,kk) < 0,kk)))))
                pondage_V(hh,kk) = 0;
            else
                pondage_V(hh,kk) = (max(DeltaV(hh,DeltaV(hh,:,kk) > 0,kk)) + abs(min(DeltaV(hh,DeltaV(hh,:,kk) < 0,kk)))).*1e-6;
            end
            % add reservoir hydro power plant
            % subtract run-of-the-river contribution
            P_load = P_load - P1(hh,:,kk);
            % redefine peak hours (F1 + F2)
            clear peak_hours
            peak_hours(1,:) = find(idx_tou(:,hh) ~= 3);
            % actual power delivered in a working day of January
            % N.B. P_load is w/o P1
            P2(hh,peak_hours,kk) = P_load(peak_hours) - P_bl;
            used_volume(hh,kk) = sum(compute_V(P2(hh,:,kk),H));
            avail_vol_2 = avail_vol_2 - used_volume(hh,kk).*no_days(hh,kk);
            % restore P_load at original value
            P_load = pu_prof(:,hh)'.*Pn;
            Pg(hh,:,kk) = P_load - P1(hh,:,kk) - P2(hh,:,kk);
            if FIGURES
                figure(4*(kk-1)+2);
                subplot(3,1,hh)
                hold on
                p(1) = stairs([0:24],horzcat(P_load(1),P_load),'LineWidth',2);
                p(2) = stairs([0:24],horzcat(P1(hh,1,kk),P1(hh,:,kk)),'LineWidth',2);
                p(3) = stairs([0:24],horzcat(P2(hh,1,kk),P2(hh,:,kk)),'LineWidth',2);
                p(4) = stairs([0:24],horzcat(Pg(hh,1,kk),Pg(hh,:,kk)),'LineWidth',2);
                set(gca,'XLim',[0 24],'XTick',[0:4:24],'YLim',[0 180],'FontSize',11,'FontWeight','demi');
                xlabel('time [h]')
                ylabel('Power [MW]')
                legend('Load','P_1','P_2','P_g','Location','best')
                grid on
                box on
                hold off
                figure(4*(kk-1)+3);
                subplot(3,1,hh)
                hold on
                p(1) = stairs([0:24],horzcat(P_load(1),P_load),'LineWidth',2);
                overall = (P1(hh,:,kk) + P2(hh,:,kk) + Pg(hh,:,kk));
                p(2) = stairs([0:24],horzcat(overall(1),overall),'LineWidth',2);
                minusP1 = (P2(hh,:,kk) + Pg(hh,:,kk));
                p(3) = stairs([0:24],horzcat(minusP1(1),minusP1),'LineWidth',2);
                onlyPg = Pg(hh,:,kk);
                p(4) = stairs([0:24],horzcat(onlyPg(1),onlyPg),'LineWidth',2);
                set(gca,'XLim',[0 24],'XTick',[0:4:24],'YLim',[0 180],'FontSize',11,'FontWeight','demi');
                xlabel('time [h]')
                ylabel('Power [MW]')
                legend('Load','P_1','P_2','P_g','Location','SouthWest')
                title('stacked contributions')
                grid on
                box on
                hold off
                figure(4*(kk-1)+4)
                subplot(3,1,hh)
                hold on
                p(1) = stairs([0:24],horzcat(P_load(1),P_load),'LineWidth',2);
                p(2) = stairs([0:24],horzcat(overall(1),overall),'LineWidth',2);
                p(3) = stairs([0:24],horzcat(minusP1(1),minusP1),'LineWidth',2,'Color',p(3).Color);
                p(4) = stairs([0:24],horzcat(onlyPg(1),onlyPg),'LineWidth',2);
                set(gca,'XLim',[0 24],'XTick',[0:4:24],'YLim',[0 180],'FontSize',11,'FontWeight','demi')
                xlabel('time [h]')
                ylabel('Power [MW]')
                x = [p(4).XData(1),repelem(p(4).XData(2:end),2)];
                y = [repelem(p(4).YData(1:end-1),2),p(4).YData(end)];
                p(7) = fill([x,fliplr(x)],[y,zeros(size(y))], p(4).Color);
                set(p(7),'EdgeColor','none');
                set(get(get(p(7),'Annotation'),'LegendInformation'),'IconDIsplayStyle','off');
                x = [p(3).XData(1),repelem(p(3).XData(2:end),2)];
                y = [repelem(p(3).YData(1:end-1),2),p(3).YData(end)];
                p(6) = fill([x,fliplr(x)],[y,p(4).YData(end),repelem(p(4).YData((end-1):-1:1),2)], p(3).Color);
                set(p(6),'EdgeColor','none');
                set(get(get(p(6),'Annotation'),'LegendInformation'),'IconDIsplayStyle','off');
                x = [p(2).XData(1),repelem(p(2).XData(2:end),2)];
                y = [repelem(p(2).YData(1:end-1),2),p(2).YData(end)];
                p(5)= fill([x,fliplr(x)],[y,p(3).YData(end),repelem(p(3).YData((end-1):-1:1),2)], p(2).Color);
                set(p(5),'EdgeColor','none');
                set(get(get(p(5),'Annotation'),'LegendInformation'),'IconDIsplayStyle','off');
                p(8) = stairs([0:24],horzcat(P_load(1),P_load),'LineWidth',2,'Color',p(1).Color);
                set(get(get(p(8),'Annotation'),'LegendInformation'),'IconDIsplayStyle','off');
                legend('Load','P_1','P_2','P_g','Location','NorthWest')
                title('energies')
                grid on
                box on
                hold off
            end
            clear idx_ps
        end
    end
end

%%
% A working day in January
ii = 1;




% A pre-holiday in January
ii = 1;
P_load = phl_pu.*Pn;
h = figure;
hold on
stairs([0:24],vertcat(P_load(1),P_load)','LineWidth',2)
set(gca,'XLim',[0 24],'XTick',[0:4:24],'YLim',[0 180],'FontSize',11,'FontWeight','demi');
xlabel('time [h]')
ylabel('Power [MW]')
grid on
box on
hold off
% baseload
P_bl = min(P_load)
P(q1(ii),Hm)
% peak hours and max load
clear idx_ps
[~,idx_ps(1,:)] = find(diag((ph_tou == 2) & (phl_pu > .55)) == 1);
peak_start = idx_ps(1)
peak_end = idx_ps(end)
peak_hours = [peak_start:1:peak_end]
E_ph_tot = sum(phl_pu(peak_hours).*Pn-P_bl.*ones(size(peak_hours)))
% run-of-the-river hydro
P1 = zeros(1,24);
% total available energy in a day (with 20 m^3/s)
E1_av = P(q1(ii),Hm).*24
idx_p1 = peak_hours
% average power and assignment
if E1_av >E_ph_tot
    E_ph_tot./length(idx_p1)
    P1(idx_p1) = E_ph_tot./length(idx_p1);
    P1(~ismember([1:24],idx_p1)) = (E1_av - E_ph_tot)./24;
    P1(idx_p1) = P1(idx_p1) + (E1_av - E_ph_tot)./24;
else
    E1_av./length(idx_p1)
    P1(idx_p1) = E1_av./length(idx_p1);
end
% assign this profile to P1

% compute the profiles of the volumes
% volumes of the controlled inflow
ctrl_infl = cumsum(compute_V(P1,Hm)); % ctrl_infl = cumsum(compute_V(P1(hh:kk),Hm));
% volumes of the natural inflow
nat_infl = cumsum(q1(kk).*ones(1,24).*3600);
h = figure;
hold on
p(1) = plot([0:24],horzcat(0,nat_infl).*1e-6,'LineWidth',2);
p(2) = plot([0:24],horzcat(0,ctrl_infl).*1e-6,'LineWidth',2);
set(gca,'XLim',[0 24],'XTick',[0:4:24],'FontSize',11,'FontWeight','demi')
xlabel('time [h]')
ylabel('volume [Mm^3]')
legend('natural inflow','controlled inflow','Location','best')
grid on
box on
hold off
% difference between natural inflow straight line and controlled inflow
% curve
DeltaV = horzcat(0,nat_infl) - horzcat(0,ctrl_infl);
% volume of the pondage in Mm^3
pondage_V = (max(DeltaV(DeltaV > 0)) + abs(min(DeltaV(DeltaV < 0)))).*1e-6
% add reservoir hydro power plant
% subtract run-of-the-river contribution
P_load = P_load - P1;
if max(P_load < P_bl)
    P2 = zeros(1,24);
else
    sum((P_load(peak_hours)-P_bl))
    sum((P_load(peak_hours)-P_bl))./24
    % difference w.r.t. average power
    P2_ave - sum((P_load(peak_hours)-P_bl))./24
    % actual power delivered in a pre-holiday of January
    % N.B. P_load is w/o P1
    P2 = zeros(1,24);
    P2(peak_hours) = P_load(peak_hours) - P_bl;
    avail_vol_2 = avail_vol_2 - sum(compute_V(P2,H));
    avail_vol_2.*1e-6
end
% average power (in a day)

% restore P_load at original value
P_load = phl_pu.*Pn;
Pg = P_load'-P1 - P2;
h = figure;
hold on
p(1) = stairs([0:24],vertcat(P_load(1),P_load),'LineWidth',2);
p(2) = stairs([0:24],horzcat(P1(1),P1),'LineWidth',2);
p(3) = stairs([0:24],horzcat(P2(1),P2),'LineWidth',2);
p(4) = stairs([0:24],horzcat(Pg(1),Pg),'LineWidth',2);
set(gca,'XLim',[0 24],'XTick',[0:4:24],'YLim',[0 180],'FontSize',11,'FontWeight','demi');
xlabel('time [h]')
ylabel('Power [MW]')
legend('Load','P_1','P_2','P_g','Location','best')
grid on
box on
hold off
h = figure;
hold on
p(1) = stairs([0:24],vertcat(P_load(1),P_load),'LineWidth',2);
p(2) = stairs([0:24],horzcat(P1(1)+P2(1)+Pg(1),P1+P2+Pg),'LineWidth',2);
p(3) = stairs([0:24],horzcat(P2(1)+Pg(1),P2+Pg),'LineWidth',2);
p(4) = stairs([0:24],horzcat(Pg(1),Pg),'LineWidth',2);
set(gca,'XLim',[0 24],'XTick',[0:4:24],'YLim',[0 180],'FontSize',11,'FontWeight','demi');
xlabel('time [h]')
ylabel('Power [MW]')
legend('Load','P_1','P_2','P_g','Location','SouthWest')
title('stacked contributions')
grid on
box on
hold off
h = figure;
hold on
p(1) = stairs([0:24],vertcat(P_load(1),P_load),'LineWidth',2);
p(2) = stairs([0:24],horzcat(P2(1)+P1(1)+Pg(1),P1+P2+Pg),'LineWidth',2);
p(3) = stairs([0:24],horzcat(P2(1)+Pg(1),P2+Pg),'LineWidth',2,'Color',p(3).Color);
p(4) = stairs([0:24],horzcat(Pg(1),Pg),'LineWidth',2);
set(gca,'XLim',[0 24],'XTick',[0:4:24],'YLim',[0 180],'FontSize',11,'FontWeight','demi')
xlabel('time [h]')
ylabel('Power [MW]')


x = [p(4).XData(1),repelem(p(4).XData(2:end),2)];
y = [repelem(p(4).YData(1:end-1),2),p(4).YData(end)];
p(7) = fill([x,fliplr(x)],[y,zeros(size(y))], p(4).Color);
set(p(7),'EdgeColor','none');
set(get(get(p(7),'Annotation'),'LegendInformation'),'IconDIsplayStyle','off');
x = [p(3).XData(1),repelem(p(3).XData(2:end),2)];
y = [repelem(p(3).YData(1:end-1),2),p(3).YData(end)];
x1 = [p(3).XData(1),p(3).XData,p(3).XData(end)];
y1 = [p(3).YData(1),p(3).YData,p(3).YData(end)];
p(6) = fill([x,fliplr(x)],[y,p(4).YData(end),repelem(p(4).YData((end-1):-1:1),2)], p(3).Color);
set(p(6),'EdgeColor','none');
set(get(get(p(6),'Annotation'),'LegendInformation'),'IconDIsplayStyle','off');
x = [p(2).XData(1),repelem(p(2).XData(2:end),2)];
y = [repelem(p(2).YData(1:end-1),2),p(2).YData(end)];
p(5)= fill([x,fliplr(x)],[y,p(3).YData(end),repelem(p(3).YData((end-1):-1:1),2)], p(2).Color);
set(p(5),'EdgeColor','none');
set(get(get(p(5),'Annotation'),'LegendInformation'),'IconDIsplayStyle','off');
p(8) = stairs([0:24],vertcat(P_load(1),P_load),'LineWidth',2,'Color',p(1).Color);
set(get(get(p(8),'Annotation'),'LegendInformation'),'IconDIsplayStyle','off');
legend('Load','P_1','P_2','P_g','Location','NorthWest')
title('energies')
grid on
box on
hold off

%A holiday in January
ii = 5;
P_load = hl_pu.*Pn;
h = figure;
hold on
stairs([0:24],vertcat(P_load(1),P_load),'LineWidth',2)
set(gca,'XLim',[0 24],'XTick',[0:4:24],'YLim',[0 180],'FontSize',11,'FontWeight','demi');
xlabel('time [h]')
ylabel('Power [MW]')
grid on
box on
hold off
% baseload
P(q1(ii),Hm)
% peak hours and max load
clear idx_ps
idx_ps(1,:) = [1:24];
peak_start = idx_ps(1)
peak_end = idx_ps(end)
peak_hours = [peak_start:1:peak_end];
E_ph_tot = sum(hl_pu(peak_hours).*Pn)
% run-of-the-river hydro
P1 = zeros(1,24);
% total available energy in a day (with 20 m^3/s)
E1_av = P(q1(ii),Hm).*24
idx_p1 = peak_hours;
% average power
E1_av./length(idx_p1)
% assign this profile to P1
P1(idx_p1) = E1_av./length(idx_p1);
% the baseload is the difference between load and this profile
P_bl = P_load(1)-P1(1)
% compute the profiles of the volumes
% volumes of the controlled inflow
ctrl_infl = cumsum(compute_V(P1,Hm));
% volumes of the natural inflow
nat_infl = cumsum(q1(ii).*ones(1,24).*3600);
h = figure;
hold on
p(1) = plot([0:24],horzcat(0,nat_infl).*1e-6,'LineWidth',2);
p(2) = plot([0:24],horzcat(0,ctrl_infl).*1e-6,'LineWidth',2);
set(gca,'XLim',[0 24],'XTick',[0:4:24],'FontSize',11,'FontWeight','demi')
xlabel('time [h]')
ylabel('volume [Mm^3]')
legend('natural inflow','controlled inflow','Location','best')
grid on
box on
hold off
% difference between natural inflow straight line and controlled inflow
% curve
DeltaV = horzcat(0,nat_infl) - horzcat(0,ctrl_infl);
% volume of the pondage in Mm^3
if isempty((max(DeltaV(DeltaV > 0)) + abs(min(DeltaV(DeltaV < 0)))))
    pondage_V = 0
else
    pondage_V = (max(DeltaV(DeltaV > 0)) + abs(min(DeltaV(DeltaV < 0)))).*1e-6
end
% add reservoir hydro power plant
% P2 is not used
P2 = zeros(1,24);
h = figure;
hold on
p(1) = stairs([0:24],vertcat(P_load(1),P_load),'LineWidth',2);
p(2) = stairs([0:24],horzcat(P1(1),P1),'LineWidth',2);
p(3) = stairs([0:24],horzcat(P2(1),P2),'LineWidth',2);
set(gca,'XLim',[0 24],'XTick',[0:4:24],'YLim',[0 180],'FontSize',11,'FontWeight','demi');
xlabel('time [h]')
ylabel('Power [MW]')
legend('Load','P_1','P_2','Location','best')
grid on
box on
hold off
h = figure;
hold on
p(1) = stairs([0:24],vertcat(P_load(1),P_load),'LineWidth',2);
p(2) = stairs([0:24],horzcat(P1(1)+P2(1),P1+P2)+P_bl,'LineWidth',2);
p(3) = stairs([0:24],horzcat(P2(1),P2)+P_bl,'LineWidth',2);
p(4) = stairs([0:24],P_bl*ones(1,25),'LineWidth',2);
set(gca,'XLim',[0 24],'XTick',[0:4:24],'YLim',[0 180],'FontSize',11,'FontWeight','demi');
xlabel('time [h]')
ylabel('Power [MW]')
legend('Load','P_1','P_2','P_g','Location','SouthWest')
title('stacked contributions')
grid on
box on
hold off
h = figure;
hold on
p(1) = stairs([0:24],vertcat(P_load(1),P_load),'LineWidth',2);
p(2) = stairs([0:24],horzcat(P2(1)+P1(1),P1+P2)+P_bl,'LineWidth',2);
p(3) = stairs([0:24],horzcat(P2(1),P2)+P_bl,'LineWidth',2,'Color',p(3).Color);
p(4) = stairs([0:24],P_bl.*ones(1,25),'LineWidth',2);
set(gca,'XLim',[0 24],'XTick',[0:4:24],'YLim',[0 180],'FontSize',11,'FontWeight','demi')
xlabel('time [h]')
ylabel('Power [MW]')
x = [p(4).XData(1),repelem(p(4).XData(2:end),2)];
y = [repelem(p(4).YData(1:end-1),2),p(4).YData(end)];
p(7) = fill([x,fliplr(x)],[y,zeros(size(y))], p(4).Color);
set(p(7),'EdgeColor','none');
set(get(get(p(7),'Annotation'),'LegendInformation'),'IconDIsplayStyle','off');
x = [p(3).XData(1),repelem(p(3).XData(2:end),2)];
y = [repelem(p(3).YData(1:end-1),2),p(3).YData(end)];
p(6) = fill([x,fliplr(x)],[y,p(4).YData(end),repelem(p(4).YData((end-1):-1:1),2)], p(3).Color);
set(p(6),'EdgeColor','none');
set(get(get(p(6),'Annotation'),'LegendInformation'),'IconDIsplayStyle','off');
x = [p(2).XData(1),repelem(p(2).XData(2:end),2)];
y = [repelem(p(2).YData(1:end-1),2),p(2).YData(end)];
p(5)= fill([x,fliplr(x)],[y,p(3).YData(end),repelem(p(3).YData((end-1):-1:1),2)], p(2).Color);
set(p(5),'EdgeColor','none');
set(get(get(p(5),'Annotation'),'LegendInformation'),'IconDIsplayStyle','off');
p(8) = stairs([0:24],vertcat(P_load(1),P_load),'LineWidth',2,'Color',p(1).Color);
set(get(get(p(8),'Annotation'),'LegendInformation'),'IconDIsplayStyle','off');
legend('Load','P_1','P_2','P_g','Location','NorthWest')
title('energies')
grid on
box on
hold off
