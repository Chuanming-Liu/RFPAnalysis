function plot_v(model)

lanz=length(model.vp); % number of layers
vmod=zeros(lanz,6);
vmod(:,1)=model.z;
vmod(:,2)=model.vp;
vmod(:,3)=model.vs;
vmod(:,4)=model.rho;
vmod(:,5)=model.qp;
vmod(:,6)=model.qs;
dz=[model.conr model.moho model.d410 model.d520 model.d660 model.cmb model.icb model.dz];
[dz,sorter]=sort(dz);
dname=strvcat('Conrad','Moho','Olivine ab','Olivine bg','Olivine Perovskite','CMB','ICB',model.dname); % discontinuity names
dname=dname(sorter,:);
modelname=model.name;
%%% remove non-existing standard discontinuities
indies=find(~isnan(dz));
dz=dz(indies);
dname=dname(indies,:);


%%% add empty layer to define planets surface
vmod=[0 -1 -1 -1 -1 -1; vmod];
dz=[0 dz];
dname=strvcat('surface',dname);

%%%%% some important values
anz=length(dz); % number of discontinuities
lanz=lanz+1; % number of layers - add one for the surface!
rplanet=model.rp; % radius of planet
vsmax=max(vmod(:,3)); % max S velocity
vpmax=max(vmod(:,2)); % max P velocity

plot(vmod(:,2),vmod(:,1),'b--');
hold on
plot(vmod(:,3),vmod(:,1),'b--');
hold on
%% mark discontinuities
set(gca,'YTick',[dz rplanet]);
set(gca,'YGrid','on');
%set(gca,'XGrid','on');
%% axes and decoration
axis ij
mkaxis('y',[min(vmod(:,1)) max(vmod(:,1))]);
mkaxis('x',[0 15]);
xlabel('velocity [km/s]');
ylabel('depth [km]');
%title(['velocity model: ' modelname]);