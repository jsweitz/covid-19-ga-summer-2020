clf
clear all
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'fignearterm_0328';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);
set(gcf,'position', [514 145 687 799]);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');

% main data goes here
% loglog(,, '');
% COVID-19 Age Structured Model
% Structure has 2 layers
% Layer 1 - Free to Move
% Layer 2 - Hospitals

% Reset

% Population
agepars.meanage=5:10:95;
agepars.highage=[9:10:99];  % Age groups
agepars.lowage=[0:10:90];  % Age groups
% From 2018 Census
population.N=10666108;
population.agefrac = [0.126	0.137	0.139	0.132	0.130	0.129	0.104	0.061	0.036	0.007];
population.agefrac=population.agefrac/sum(population.agefrac); % Must sum to 1
population.meanage = sum(agepars.meanage.*population.agefrac);

% Parameters
% Parameters
pars.gamma_e=1/4;   % Transition to infectiousness
pars.gamma_a=1/6;   % Resolution rate for asymptomatic
pars.gamma_s=1/6;  % Resolution rate for symptomatic
pars.gamma_h=1/10;  % Resolution rate in hospitals
pars.beta_a=4/10;   % Transmission for asymptomatic
pars.beta_s=8/10;      % Transmission for symptomatic
% pars.p=0.5;         % Fraction asymptomatic
% pars.p=[0.99 0.99 0.95 0.9 0.8 0.7 0.6 0.5 0.5 0.5];         % Fraction asymptomatic
pars.p=[0.95 0.95 0.90 0.8 0.7 0.6 0.4 0.2 0.2 0.2];         % Fraction asymptomatic
pars.overall_p=sum(pars.p.*population.agefrac);
pars.Itrigger = 45000/population.N; % Trigger at 5000 total cases, irrespective of type


% Epi parameters
pars.Ra=pars.beta_a/pars.gamma_a;
pars.Rs=pars.beta_s/pars.gamma_s;
pars.R0=pars.p*pars.Ra+(1-pars.p)*pars.Rs;
pars.R0=sum(pars.p.*population.agefrac*pars.Ra+(1-pars.p).*population.agefrac*pars.Rs);
pars_bau=pars;

% Age-stratification
agepars.meanage=5:10:95;
agepars.highage=[9:10:99];  % Age groups
agepars.lowage=[0:10:90];  % Age groups
agepars.hosp_frac=[0.1 0.3 1.2 3.2 4.9 10.2 16.6 24.3 27.3 27.3]/100;
agepars.hosp_crit=[5 5 5 5 6.3 12.2 27.4 43.2 70.9 70.9]/100;
agepars.crit_die= 0.5*ones(size(agepars.meanage));
agepars.num_ages = length(agepars.meanage);
N=agepars.num_ages;
agepars.S_ids=1:N;
agepars.E_ids=(N+1):2*N;
agepars.Ia_ids=(2*N+1):3*N;
agepars.Is_ids=(3*N+1):4*N;
agepars.Ihsub_ids=(4*N+1):5*N;
agepars.Ihcri_ids=(5*N+1):6*N;
agepars.R_ids=(6*N+1):7*N;
agepars.D_ids=(7*N+1):8*N;
agepars.Hcum_ids=(8*N+1):9*N;
agepars.IFR = sum((1-pars.p).*population.agefrac.*agepars.hosp_frac.*agepars.hosp_crit.*agepars.crit_die);

% Init the population - baseline
% Open plus hospitals
% SEIaIS (open) and then I_ha I_hs and then R (open) and D (cumulative) age stratified
tmpzeros = zeros(size(agepars.meanage));
outbreak.y0=[population.agefrac tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros];
% Initiate an outbreak with 500 symptomatic current caseas and 7500 asymptomatic cases
% effective 8000 total and 25 deaths (based on GA estimates)
% Initiate an outbreak
outbreak.y0=population.N*outbreak.y0;
outbreak.y0(3)=outbreak.y0(3)-1;
outbreak.y0(13)=1;
outbreak.y0=outbreak.y0/population.N;
outbreak.pTime=365;
outbreak.pNear=30;
outbreak.pshift=0;

% Sims - Get to Crossing
pars.alpha=0;  % Shielding
pars.beta_a=4/10;   % Transmission for asymptomatic
pars.beta_s=8/10;      % Transmission for symptomatic
t0_opt = 67 + 7;  % From parameter fitting;
te=t0_opt;
% Trigger approach
%opts=odeset('reltol',1e-8,'maxstep',0.1,'events',@intervene_trigger);
%[tpre,ypre,te,ye,ie]=ode45(@covid_model_ga,[0:1:outbreak.pTime], outbreak.y0,opts,pars,agepars);
% Parameter fit approach
opts=odeset('reltol',1e-8,'maxstep',0.1);
[tpre,ypre]=ode45(@covid_model_ga,[0:1:t0_opt], outbreak.y0,opts,pars,agepars);
% Baseline case
[tpost,ypost]=ode45(@covid_model_ga,[0:1:outbreak.pNear],ypre(end,:),opts,pars,agepars);
tbau=[tpre(1:end-1)-te; tpost];
ybau=[ypre(1:end-1,:); ypost];

% Stats
statsbau.R=ybau(:,agepars.R_ids);
statsbau.D=ybau(:,agepars.D_ids);
statsbau.Htot=ybau(:,agepars.Ihsub_ids)+ybau(:,agepars.Ihcri_ids);
statsbau.Hacu=ybau(:,agepars.Ihcri_ids);
statsbau.Dday_age=statsbau.D(2:end,:)-statsbau.D(1:end-1,:);
statsbau.Dday=sum(statsbau.Dday_age');
statsbau.Dcum=sum(statsbau.D');
statsbau.Hacu_day=sum(statsbau.Hacu');
statsbau.fracI=1-sum(ybau(:,agepars.S_ids)');
statsbau.Is=sum(ybau(:,agepars.Is_ids)');
statsbau.Iday=statsbau.Is(2:end)-statsbau.Is(1:end-1);
statsbau.Isevere=sum(ybau(:,3*N+1:6*N)');
tmpi=find(tbau<=0);
Htot=sum(statsbau.Htot(tmpi,:)');
Dcum=statsbau.Dcum(tmpi);

% Sims - 25%
% Irrelevant, no shielding
pars.beta_a=4*0.75/10;   % Transmission for asymptomatic
pars.beta_s=8*0.75/10;      % Transmission for symptomatic
pars.Ra=pars.beta_a/pars.gamma_a;
pars.Rs=pars.beta_s/pars.gamma_s;
pars.R0=pars.p*pars.Ra+(1-pars.p)*pars.Rs;
pars.R0=sum(pars.p.*population.agefrac*pars.Ra+(1-pars.p).*population.agefrac*pars.Rs);
pars_base=pars;
opts=odeset('reltol',1e-8,'maxstep',0.1);
[tpost,ypost]=ode45(@covid_model_ga,[0:1:outbreak.pNear],ypre(end,:),opts,pars,agepars);
tb=[tpre(1:end-1)-te; tpost];
yb=[ypre(1:end-1,:); ypost];

% Stats
statsb.R=yb(:,agepars.R_ids);
statsb.D=yb(:,agepars.D_ids);
statsb.Htot=yb(:,agepars.Ihsub_ids)+yb(:,agepars.Ihcri_ids);
statsb.Hacu=yb(:,agepars.Ihcri_ids);
statsb.Dday_age=statsb.D(2:end,:)-statsb.D(1:end-1,:);
statsb.Dday=sum(statsb.Dday_age');
statsb.Dcum=sum(statsb.D');
statsb.Hacu_day=sum(statsb.Hacu');
statsb.fracI=1-sum(yb(:,agepars.S_ids)');
statsb.Is=sum(yb(:,agepars.Is_ids)');
statsb.Iday=statsb.Is(2:end)-statsb.Is(1:end-1);

% Sims - Medium, Some Control
pars.alpha=0;  % Shielding
pars.beta_a=4*0.5/10;   % Transmission for asymptomatic
pars.beta_s=8*0.5/10;      % Transmission for symptomatic
pars.Ra=pars.beta_a/pars.gamma_a;
pars.Rs=pars.beta_s/pars.gamma_s;
pars.R0=pars.p*pars.Ra+(1-pars.p)*pars.Rs;
pars.R0=sum(pars.p.*population.agefrac*pars.Ra+(1-pars.p).*population.agefrac*pars.Rs);
pars_med=pars;
[tpost,ypost]=ode45(@covid_model_ga,[0:1:outbreak.pNear],ypre(end,:),opts,pars,agepars);
t=[tpre(1:end-1)-te; tpost];
y=[ypre(1:end-1,:); ypost];

% Stats - 40%
stats.R=y(:,agepars.R_ids);
stats.D=y(:,agepars.D_ids);
stats.Htot=y(:,agepars.Ihsub_ids)+y(:,agepars.Ihcri_ids);
stats.Hacu=y(:,agepars.Ihcri_ids);
stats.Dday_age=stats.D(2:end,:)-stats.D(1:end-1,:);
stats.Dday=sum(stats.Dday_age');
stats.Hacu_day=sum(stats.Hacu');
stats.fracI=1-sum(y(:,agepars.S_ids)');
stats.Is=sum(y(:,agepars.Is_ids)');
stats.Iday=stats.Is(2:end)-stats.Is(1:end-1);
stats.Dcum=sum(stats.D');


% Sims - Low, Severe Control
pars.beta_a=4*0.25/10;   % Transmission for asymptomatic
pars.beta_s=8*0.25/10;      % Transmission for symptomatic
pars.Ra=pars.beta_a/pars.gamma_a;
pars.Rs=pars.beta_s/pars.gamma_s;
pars.R0=pars.p*pars.Ra+(1-pars.p)*pars.Rs;
pars.R0=sum(pars.p.*population.agefrac*pars.Ra+(1-pars.p).*population.agefrac*pars.Rs);
pars_low=pars;
[tpost,ypost]=ode45(@covid_model_ga,[0:1:outbreak.pNear],ypre(end,:),opts,pars,agepars);
th=[tpre(1:end-1)-te; tpost];
yh=[ypre(1:end-1,:); ypost];

% Stats
statsh.R=yh(:,agepars.R_ids);
statsh.D=yh(:,agepars.D_ids);
statsh.Htot=yh(:,agepars.Ihsub_ids)+yh(:,agepars.Ihcri_ids);
statsh.Hacu=yh(:,agepars.Ihcri_ids);
statsh.Dday_age=statsh.D(2:end,:)-statsh.D(1:end-1,:);
statsh.Dday=sum(statsh.Dday_age');
statsh.Dcum=sum(statsh.D');
statsh.Hacu_day=sum(statsh.Hacu');
statsh.fracI=1-sum(yh(:,agepars.S_ids)');
statsh.Is=sum(yh(:,agepars.Is_ids)');
statsh.Iday=statsh.Is(2:end)-statsh.Is(1:end-1);

%semilogy(t,y);
%legend('S','E','I1','I2','R','D');

subplot(3,1,1);
curdate=datetime('now')-1;
tickdates = curdate+[-10:2:outbreak.pNear];
formatOut = 'mm/dd';
xticklabels=datestr(tickdates,formatOut);
tmpi=find(t>=0);
tmph=plot(curdate+tbau(tmpi),statsbau.Dcum(tmpi)*population.N,'k-');
set(tmph,'linewidth',3);
hold on
tmph=plot(curdate+t(tmpi),statsb.Dcum(tmpi)*population.N,'b--');
set(tmph,'linewidth',3);
tmph=plot(curdate+t(tmpi),stats.Dcum(tmpi)*population.N,'g:');
set(tmph,'linewidth',3);
tmph=plot(curdate+t(tmpi),statsh.Dcum(tmpi)*population.N,'ko');
set(tmph,'markersize',8,'markerfacecolor','g');
tmpi=find(t<=0);
tmph=plot(curdate+t(tmpi),statsh.Dcum(tmpi)*population.N,'k-');
set(tmph,'linewidth',3,'color',[0.5 0.5 0.5]);
set(gca,'fontsize',16);
set(gca,'xtick',tickdates,'xticklabelrotation',30,'fontsize',12);
set(gca,'xticklabels',xticklabels);
xlim([tickdates(1) tickdates(end)]);
xlabel('Date','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('Cumulative deaths','fontsize',18,'verticalalignment','bottom','interpreter','latex');
title({'COVID-19 Near-Term Epidemic, Georgia Assessment, March 28, 2020';'Updated 1 Month Projection, Age-Stratified Risk, No Prior Distancing Effects'},'fontsize',18,'interpreter','latex')
ylim([0 4000]);
set(gca,'ytick',[0:500:4000]);
%set(gca,'yticklabels',{'0';'5,000';'10,000';'15,000';'20,000';'25,000';'30,000'});
tmpl=legend('Baseline','25\% Contact Reduction','50\% Contact Reduction','75\% Contact Reduction');
set(tmpl,'interpreter','latex','location','northwest','fontsize',14);
legend('boxoff');

subplot(3,1,2);
tmpi=find(t>=0);
tmph=plot(curdate+t(tmpi),statsbau.Hacu_day(tmpi)*population.N,'k-');
set(tmph,'linewidth',3);
hold on
tmph=plot(curdate+t(tmpi),statsb.Hacu_day(tmpi)*population.N,'b--');
set(tmph,'linewidth',3);
tmph=plot(curdate+t(tmpi),stats.Hacu_day(tmpi)*population.N,'g:');
set(tmph,'linewidth',3);
tmph=plot(curdate+t(tmpi),statsh.Hacu_day(tmpi)*population.N,'ko');
set(tmph,'markersize',8,'markerfacecolor','g');
tmpi=find(t<=0);
tmph=plot(curdate+t(tmpi),statsh.Hacu_day(tmpi)*population.N,'k-');
set(tmph,'linewidth',3,'color',[0.5 0.5 0.5]);
tmph=plot([curdate-te curdate+outbreak.pTime],[3378 3378]/3,'r-');
set(tmph,'linewidth',2,'color',[0.85 0 0]);
tmpt=text(curdate-8,1350,{'ICU Surge Capacity w/33% Availability'});
set(gca,'fontsize',16);
tmph=plot([curdate+10 curdate+outbreak.pTime],[3378 3378],'r-');
set(tmph,'linewidth',2,'color',[0.85 0 0]);
tmpt=text(curdate+20,3650,{'ICU Max Surge'});
set(gca,'fontsize',16);
set(tmpt,'fontsize',14,'interpreter','latex');
set(gca,'xtick',tickdates,'xticklabelrotation',30,'fontsize',12);
set(gca,'xticklabels',xticklabels);
xlim([tickdates(1) tickdates(end)]);
xlabel('Date','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('ICU beds needed','fontsize',18,'verticalalignment','bottom','interpreter','latex');
% title('','fontsize',24)
ylim([0 4000]);
tmpl=legend('Baseline','25\% Contact Reduction','50\% Contact Reduction','75\% Contact Reduction');
set(tmpl,'interpreter','latex','location','northwest','fontsize',14);
legend('boxoff');
subplot(3,1,3);
tmpi=find(t>=0);
tmph=plot(curdate+tbau(tmpi),statsbau.fracI(tmpi)*(1-pars_bau.overall_p)*population.N,'k-');
statsbau.fracI=1-sum(ybau(:,agepars.S_ids)');
set(tmph,'linewidth',3);
hold on
tmph=plot(curdate+tb(tmpi),statsb.fracI(tmpi)*(1-pars_base.overall_p)*population.N,'b--');
set(tmph,'linewidth',3);
tmph=plot(curdate+t(tmpi),stats.fracI(tmpi)*(1-pars_med.overall_p)*population.N,'g:');
set(tmph,'linewidth',3);
tmph=plot(curdate+th(tmpi),statsh.fracI(tmpi)*(1-pars_low.overall_p)*population.N,'ko');
set(tmph,'markersize',8,'markerfacecolor','g');
tmpi=find(t<=0);
tmph=plot(curdate+tbau(tmpi),statsbau.fracI(tmpi)*(1-pars_bau.overall_p)*population.N,'k-');
set(tmph,'linewidth',3,'color',[0.5 0.5 0.5]);
set(gca,'fontsize',16);
set(gca,'xtick',tickdates,'xticklabelrotation',30,'fontsize',12);
set(gca,'xticklabels',xticklabels);
set(gca,'ytick',[0:20000:100000]);
set(gca,'yticklabels',{'0';'20,000';'40,000';'60,000';'80,000';'100,000'});
xlim([tickdates(1) tickdates(end)]);
xlabel('Date','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel({'Cumulative';'symptomatic cases'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
% title('','fontsize',24)
ylim([0 10^5]);
tmpl=legend('Baseline','25\% Contact Reduction','50\% Contact Reduction','75\% Contact Reduction');
set(tmpl,'interpreter','latex','location','northwest','fontsize',14);
legend('boxoff');
% title('','fontsize',24)


%
%
% Some helpful plot commands
% tmph=plot(x,y,'ko');
% set(tmph,'markersize',10,'markerfacecolor,'k');
% tmph=plot(x,y,'k-');
% set(tmph,'linewidth',2);

% for use with layered plots
% set(gca,'box','off')

% adjust limits
% tmpv = axis;
% axis([]);
% ylim([]);
% xlim([]);

% change axis line width (default is 0.5)
% set(tmpa1,'linewidth',2)

% fix up tickmarks
% set(gca,'xtick',[1 100 10^4])
% set(gca,'ytick',[1 100 10^4])

% creation of postscript for papers
% psprint(tmpxfigfilename);

% the following will usually not be printed 
% in good copy for papers
% (except for legend without labels)

% legend
% tmplh = legend('stuff',...);
% tmplh = legend('','','');
% remove box
% set(tmplh,'visible','off')
% legend('boxoff');

% title('','fontsize',24)
% 'horizontalalignment','left');

% for writing over the top
% coordinates are normalized again to (0,1.0)
tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');
% first two points are normalized x, y positions
% text(,,'','Fontsize',14);

%tmpt=text(-0.2,-0.25,{'Calculation - J.S.Weitz - jsweitz@gatech.edu - 3/28/20 - Nonlinear Population Dynamics Using Age-Dependent Risk';'License: Creative Commons BY-SA 4.0, i.e., Share, Adapt, Attribute - https://creativecommons.org/licenses/by/4.0/';'Thanks to C.~Andris, K.~Carden, J~Dushoff, and Weitz group members, code https://github.com/jsweitz/covid-19-ga-summer-2020'});
%set(tmpt,'fontsize',10,'interpreter','latex');

% automatic creation of postscript
% without name/date
psprintc(tmpfilenoname);
psprint(tmpfilebwname);

tmpt = pwd;
tmpnamememo = sprintf('[source=%s/%s.ps]',tmpt,tmpprintname);
text(1.05,.05,tmpnamememo,'Fontsize',6,'rotation',90);
datenamer(1.1,.05,90);
% datename(.5,.05);
% datename2(.5,.05); % 2 rows

% automatic creation of postscript
psprintc(tmpfilename);

% set following on if zooming of 
% plots is required
% may need to get legend up as well
%axes(tmpa1)
%axes(tmplh)
clear tmp*
