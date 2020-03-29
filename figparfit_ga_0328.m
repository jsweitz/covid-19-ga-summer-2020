clf;
clear all
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'figparfit_ga_0328';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[427 539 1076 394]);

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
pars.Itrigger = 500000/population.N; % Trigger at 5000 total cases, irrespective of type


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
opts=odeset('reltol',1e-8,'maxstep',0.1,'events',@intervene_trigger);
[tpre,ypre,te,ye,ie]=ode45(@covid_model_ga,[0:1:outbreak.pTime], outbreak.y0,opts,pars,agepars);
% Get stats of sims
stats.Hcum = sum(ypre(:,agepars.Hcum_ids)');
stats.Dcum = sum(ypre(:,agepars.D_ids)');

% Timing
curdate=datetime('now');
tickdates = curdate+[-te:10:10];
formatOut = 'mm/dd';

% Data - Through 3/28 - last 7 days
ga_cases  = [552 620 800 1097 1387 1643 2198];
ga_hosp =  [186 240 361 438 509 607 660];
ga_death = [25 25 38 47 56 65 79];
ga_date = [curdate-6:curdate];


subplot(1,2,1);
tmph=semilogy(ga_date,ga_hosp,'bd');
set(tmph,'markersize',12,'markerfacecolor',[0.5 0.5 0.5],'linewidth',1);
hold on
tmph=semilogy(ga_date,ga_death,'ro');
set(tmph,'markersize',12,'markerfacecolor',[0.5 0.5 0.5],'linewidth',1);
tmph=semilogy(curdate-te+tpre,stats.Hcum*population.N,'b-');
set(tmph,'linewidth',2);
tmph=semilogy(curdate-te+tpre,stats.Dcum*population.N,'r-');
set(tmph,'linewidth',2);
ylim([10^0 10^3]);

xticklabels=datestr(tickdates,formatOut);
set(gca,'fontsize',16);
set(gca,'xtick',tickdates,'xticklabelrotation',30,'fontsize',12);
set(gca,'xticklabels',xticklabels);
xlim([tickdates(1) tickdates(end)]);
xlabel('Date','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('Cumulative counts','fontsize',18,'verticalalignment','bottom','interpreter','latex');
title('Free-running epidemic fits','fontsize',18,'interpreter','latex')
ylim([10^0 10^4]);

tmpl=legend('Hospitalization - Data','Deaths - Data','Hospitalization - Model','Deaths - Model');
set(tmpl,'interpreter','latex','location','northwest','fontsize',14);
legend('boxoff');

subplot(1,2,2);
% T-equivalent - find the equivalent time given the
% observed case data
clear sse
pars.fit_days=7;
tstart_max = floor(te)-10;
tstart_min = 20;
for t0=tstart_min:tstart_max,
   model_hosp = stats.Hcum(t0:t0+pars.fit_days-1)*population.N;
   model_death = stats.Dcum(t0:t0+pars.fit_days-1)*population.N;
   sse(t0) = sqrt(sum(log(model_hosp)-log(ga_hosp)).^2+sum(log(model_death)-log(ga_death)).^2);
end
tmpi=find(sse>0);
[ssemin tminid]=min(sse(tmpi));
t0_opt = tmpi(tminid);

% Plot the time-shifted results
tmph=semilogy(ga_date,ga_hosp,'bd');
set(tmph,'markersize',12,'markerfacecolor',[0.5 0.5 0.5],'linewidth',1);
hold on
tmph=semilogy(ga_date,ga_death,'ro');
set(tmph,'markersize',12,'markerfacecolor',[0.5 0.5 0.5],'linewidth',1);
tmph=semilogy(curdate-te+tpre+(te-t0_opt-7),stats.Hcum*population.N,'b-');
set(tmph,'linewidth',2);
tmph=semilogy(curdate-te+tpre+(te-t0_opt-7),stats.Dcum*population.N,'r-');
set(tmph,'linewidth',2);
ylim([10^0 10^3]);
xticklabels=datestr(tickdates,formatOut);
set(gca,'fontsize',16);
set(gca,'xtick',tickdates,'xticklabelrotation',30,'fontsize',12);
set(gca,'xticklabels',xticklabels);
xlim([tickdates(1) tickdates(end)]);
xlabel('Date','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('Cumulative counts','fontsize',18,'verticalalignment','bottom','interpreter','latex');
title('Time-shifted epidemic fits','fontsize',18,'interpreter','latex')
ylim([10^0 10^4]);
tmpl=legend('Hospitalization - Data','Deaths - Data','Hospitalization - Model','Deaths - Model');
set(tmpl,'interpreter','latex','location','northwest','fontsize',14);
legend('boxoff');


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
