% A transient multi-phase model with solute diffusion
% Explicit scheme
% Unknowns:
% pc, pb_N, pf_1, vc, vb, vf, thetan, vn, thetac, cc, cb, cf

%    -------------------------------------------------------------
%      back channel (b)    |    cell (c)    |    front channel (f)
%    -------------------------------------------------------------
%  x = 0                 x = xb           x = xf               x = L

% vb = dxb/dt
% vf = dxf/dt
% pb_N: pressure in the back channel at x = xb
% pf_1: pressure in the front channel at x = xf

% Length scale: micron
% To match the same scale of pressure in Pa, the mass scale is mg.
% i.e., 1 kg/m/s^2 = 1 mg/micron/s^2
% Time scale: s

% Allow actin polymerization happens on either end of the cell

% Contact: Yizeng Li (liyizeng52@hotmail.com)


clear
clc

%% Control parameters
L = 200;                % (micron) channel length
b = 3. ;                % (micron) channel width (the smaller dimension)
w = 10;                 % (micron) channel width (the larger dimension)
a = 1 - 192*b/pi^5/w;   % a correction factor for calculing the hydraulic resistance

R = 8.31451;            % (J/mol K) Ideal gas constant
T = 310;                % (K) absolute temperature

dt = 0.01;             % (s) time grid
T_total = (4/3)*3600;       % (s) Total time of simultion
dx_ref = 4;             % (micron) reference grid size 
dt_plot = 4*60;         % (s) Time interval for plotting
T_rev = (2/3)*3600;       % (s) Time for reversing cell migration
T_rev_scale = 50;      % (s) Time scaling factor for repolarization of actin and ion channels

Nt = round(T_total/dt) + 1; % total number of time steps. 
Time = linspace(0,T_total,Nt);
Nt_rev = round(T_rev/dt) + 1;
Nt_rev_scale = round(T_rev_scale/dt) + 1;

cb0 = 340;              % (mM) solute concentration at infinity (back channel)
cf0 = 340;              % (mM) solute concentration at infinity (front channel)
pb0 = 0;                % (Pa) hydrostatic pressure at infinity (back channel)
pf0 = 0;                % (Pa) hydrostatic pressure at infinity (front channel)
fextf = 0d2;            % (Pa) external force per unit area at the front of the cell
fextb = 0d2;            % (Pa) external force per unit area at the back of the cell
Thetac = 0.2;          % (mM) reference value of G-actin
Thetan = 0.2;          % (mM) reference value of F-actin
Thetat = Thetac + Thetan;   % (mM) reference total G- and F-actin, \int (thetan + thetac)dx/Lc

%% Cell dynamics parameters
ksigman = 8d3;          % (Pa /mM) Coefficient of passive actin network pressure
ksigmaa = 0*4d2;          % (Pa /mM) Coefficient of active actin network contraction
gamma = 3.0d-3;        % (1/s) constant rate of actin depolymerization
Jactinf0 = 3.5d-2;      % (micron mM/s) Jactinf = Jactinf0*thetac^f/(thetacc + thetac^f)
Jactinb0 = -0d-3;       % (micron mM/s) Constant actin (de)polymerization at the back
thetacc  = 0.2d-1;      % (mM) Critical value for actin polymerization

eta_it = 1d-3;          % (Pa s/micron^2/mM) interfacial coefficient between actin and cytosol
eta_st = 6d2;           % (Pa s/micron^2/mM) coefficient of focal adhesion between the actin and the substrate
muc = 2d-3;             % (Pa s) fluid dynamic viscosity in the cell
mub = 1d-3;             % (Pa s) fluid dynamic viscosity in the back channel
muf = 1d-3;             % (Pa s) fluid dynamic viscosity in the front channel
mum = 1d-0;             % (Pa s) membrane-cortex viscosity (tension coefficient)
kct = 5d4;              % (Pa micron) effective membrane-cortex stiffness
kad = 2d4;              % (Pa s/micron) adhesive force, Fad^b = kad*v0

Dcc = 1.d2;             % (micron^2/s) diffusion constant for solutes in the cytosol
Dcb = 1.d2;             % (micron^2/s) diffusion constant for solutes in the back channel
Dcf = 1.d2;             % (micron^2/s) diffusion constant for solutes in the front channel
Dtc = 1.d1;             % (micron^2/s) diffusion constant for theta_c

alphaf = 1.d-4;         % (micron/Pa/s) coefficient of water permeation
alphab = 1.d-4;         % (micron/Pa/s) coefficient of water permeation
gf = 5d1;               % (micron/s) passive channel coefficient at the front
gb = 5d1;               % (micron/s) passive channel coefficient at the back
Jcactiveb = 1d1;        % (mM micron/s) active flux at the back
Jcactivef = 0d1;        % (mM micron/s) active flux at the front

%% Initialization
xb = NaN(Nt,1);         % back cell position vector
xf = NaN(Nt,1);         % front cell position vector
Lb = NaN(Nt,1);         % back channel length vector
Lc = NaN(Nt,1);         % cell length vector
Lf = NaN(Nt,1);         % front chennel length vector
vc = NaN(Nt,1);         % cytosol velocity (constant in a 1D channel)
vb = NaN(Nt,1);         % back cell bounary velocity vector
vf = NaN(Nt,1);         % front cell bundary velocity vector
Jwaterb = NaN(Nt,1);
Jwaterf = NaN(Nt,1);
pb_N = NaN(Nt,1);       % pressure on the cell in the back channel
pf_1 = NaN(Nt,1);       % pressure on the cell in the front channel
ActRatio = NaN(Nt,1);   % F-actin ratio from upstream to downstream
ActInt = NaN(Nt,1);
ActCon = NaN(Nt,1);

%% Initial Conditions
Lc(1) = 50;             % (micron)
xb(1) = 50;             % (micron)
xf(1) = xb(1) + Lc(1);
Lb(1) = xb(1);
Lf(1) = L - xf(1);
vb(1) = 0;
vf(1) = 0;

nb = round(Lb(1)/dx_ref);   % Number of elements in the back channel
nc = round(Lc(1)/dx_ref);   % Number of elements in the cell
nf = round(Lf(1)/dx_ref);   % Number of elements in the front channel

dxb = Lb(1)/nb;             % Grid size in the back channel
dxc = Lc(1)/nc;             % Grid size in the cell
dxf = Lf(1)/nf;             % Grid size in the front channel

Nb = nb + 1;                % Number of nodes in the back channel
Nc = nc + 1;                % Number of nodes in the cell
Nf = nf + 1;                % Number of nodes in the front channel

vc(1) = 0;
Jwaterf(1) = 0;
Jwaterb(1) = 0;
vn = zeros(Nc,1);
thetac = Thetac*ones(Nc,1);
thetan = Thetan*ones(Nc,1);
sigma = ksigman*thetan;

pb_N(1) = pb0;
pf_1(1) = pf0;
cb = cb0*ones(Nb,1);
cc = (cb0 + cf0)/2*ones(Nc,1) - sigma/R/T;
cf = cf0*ones(Nf,1);

%% Time iteration
iplot = 0;
for it = 2:Nt
    
    f_on  = 1/(1+exp(-(it-Nt_rev)/Nt_rev_scale));
    f_off = 1/(1+exp((it-Nt_rev)/Nt_rev_scale));
    
    f_on_neg  = 1/(1+exp(-(it-Nt_rev)/Nt_rev_scale))-0.8;
    f_off_neg = 1/(1+exp((it-Nt_rev)/Nt_rev_scale))-0.2;

    f_half = 1/(1+exp((it-Nt_rev)/Nt_rev_scale))+1;
    f_third = 1/(1+exp((it-Nt_rev)/Nt_rev_scale))+0.5;
    f_fifth = 1/(1+exp((it-Nt_rev)/Nt_rev_scale))+0.25;
    f_11 = 1/(1+exp((it-Nt_rev)/Nt_rev_scale))+0.1;
    f_double = 1/(1+exp(-(it-Nt_rev)/Nt_rev_scale))+1;
    f_triple = 1/(1+exp(-(it-Nt_rev)/Nt_rev_scale))+0.5;

    pb0 = 160*f_on;

    Jactinf0 = 3.5d-2*f_off;
    Jactinb0 = 3.5d-2*f_on;
    gamma = 3.0d-3*f_third;
    Jcactivef = 1d1*f_off;
    Jcactiveb = 0d0*f_on;


    %% Water module
    % pc (Nc elements), pb(i=Nb), pf(i=1), vc, vb_water, vf_water)
    Mat = zeros(Nc+5, Nc+5);
    RHS = zeros(Nc+5, 1);
    
    npb = Nc + 1;   % index for pb(i=Nb)
    npf = Nc + 2;   % index for pf(i=1)
    nvc = Nc + 3;   % index for vc
    nvb = Nc + 4;   % index for vb
    nvf = Nc + 5;   % index for vf
    
    Mat(1, [1, 2]) = [-1 1]/dxc;
    Mat(1, nvc) = 12/a*muc/b^2 + eta_it*thetan(1);
    RHS(1) = eta_it*thetan(1)*vn(1);
    for in = 2:Nc-1
        Mat(in, [in-1, in+1]) = [-1 1]/2/dxc;
        Mat(in, nvc) = 12/a*muc/b^2 + eta_it*thetan(in);
        RHS(in) = eta_it*thetan(in)*vn(in);
    end
    Mat(Nc,   [npb, nvc]) = [1,  12/a*mub*Lb(it-1)/b^2];
    Mat(Nc+1, [npf, nvc]) = [1, -12/a*muf*Lf(it-1)/b^2];
    RHS([Nc, Nc+1]) = [pb0, pf0];
    Mat(Nc+2, [nvf, nvc, Nc, npf]) = [1, -1, alphaf, -alphaf];
    Mat(Nc+3, [nvb, nvc, 1,  npb]) = [1, -1, -alphab, alphab];
    RHS(Nc+2) =   alphaf*R*T*(cc(Nc) - cf(1));
    RHS(Nc+3) = - alphab*R*T*(cc(1)  - cb(Nb));
    Mat(Nc+4, [npf, Nc, nvf]) = [b, -b,  (2*mum + kad)];
    Mat(Nc+5, [npb, 1,  nvb]) = [b, -b, -(2*mum + kad)];
    RHS(Nc+4) = b*sigma(Nc) - 2*kct*((Lc(it-1)-Lc(1))/Lc(1)) + b*fextf;
    RHS(Nc+5) = b*sigma(1)  - 2*kct*((Lc(it-1)-Lc(1))/Lc(1)) + b*fextb;
    
    X = Mat\RHS;
    pc = X(1:Nc);
    pb_N(it) = X(npb);
    pf_1(it) = X(npf);
    vc(it) = X(nvc);
    vb(it) = X(nvb);
    vf(it) = X(nvf);
    
    % post-process
    Jwaterf(it) = -alphaf*((pc(Nc)-pf_1(it)) - R*T*(cc(Nc)-cf(1)));
    Jwaterb(it) = -alphab*((pc(1)-pb_N(it))  - R*T*(cc(1)-cb(Nb)));

    %% Location and grid update
    xb(it) = xb(it-1) + vb(it)*dt;
    xf(it) = xf(it-1) + vf(it)*dt;
    Lb(it) = xb(it);
    Lf(it) = L - xf(it);
    Lc(it) = xf(it) - xb(it);

    nb = round(Lb(it)/dx_ref);   % Number of elements in the back channel
    nc = round(Lc(it)/dx_ref);   % Number of elements in the cell
    nf = round(Lf(it)/dx_ref);   % Number of elements in the front channel

    dxb = Lb(it)/nb;             % Grid size in the back channel
    dxc = Lc(it)/nc;             % Grid size in the cell
    dxf = Lf(it)/nf;             % Grid size in the front channel

    Nb = nb + 1;                % Number of nodes in the back channel
    Nc = nc + 1;                % Number of nodes in the cell
    Nf = nf + 1;                % Number of nodes in the front channel
    
    thetan = GridUpdate_con2(thetan,xb(it-1),xf(it-1),xb(it),xf(it),dx_ref,vb(it),vf(it),dt);
    thetac = GridUpdate_con2(thetac,xb(it-1),xf(it-1),xb(it),xf(it),dx_ref,vb(it),vf(it),dt);
    cc = GridUpdate_con2(cc,xb(it-1),xf(it-1),xb(it),xf(it),dx_ref,vb(it),vf(it),dt);
    cb = GridUpdate_con2(cb,0,xb(it-1),0,xb(it),dx_ref,0,vb(it),dt);
    cf = GridUpdate_con2(cf,xf(it-1),L,xf(it),L,dx_ref,vf(it),0,dt);
    
    %% Actin module
    Jactinf = Jactinf0*thetac(Nc)/(thetacc+thetac(Nc));
    Jactinb = Jactinb0*thetac(1)/(thetacc+thetac(1));
    
    % ghost points for thetan
    thetan_ghost_b = 2/ksigman*((eta_it+eta_st)*vb(it)*dxc*thetan(1) ...
        - eta_it*vc(it)*dxc*thetan(1) + ksigman/2*thetan(2) ...
        + (eta_it+eta_st)*dxc*Jactinb);
    thetan_ghost_f = -2/ksigman*((eta_it+eta_st)*vf(it)*dxc*thetan(Nc) ...
        - eta_it*vc(it)*dxc*thetan(Nc) - ksigman/2*thetan(Nc-1) ...
        - (eta_it+eta_st)*dxc*Jactinf);
    thetan_temp = [thetan_ghost_b; thetan; thetan_ghost_f];
    % time iteration for thetan
    thetan = thetan + dt*(ksigman/(eta_it+eta_st)/dxc^2 ...
        *(thetan_temp(3:Nc+2)-2*thetan_temp(2:Nc+1)+thetan_temp(1:Nc)) ...
        - eta_it/(eta_it+eta_st)/2/dxc*vc(it)*(thetan_temp(3:Nc+2)-thetan_temp(1:Nc)) ...
        - gamma.*thetan);
    
    % ghost points for thetac
    thetac_ghost_b = thetac(2)    + 2*dxc/Dtc*(vb(it)-vc(it))*thetac(1)  - 2*dxc/Dtc*Jactinb;
    thetac_ghost_f = thetac(Nc-1) + 2*dxc/Dtc*(vc(it)-vf(it))*thetac(Nc) - 2*dxc/Dtc*Jactinf;
    thetac_temp = [thetac_ghost_b; thetac; thetac_ghost_f];
    % time iteration for thetac
    thetac = thetac + dt*(Dtc/dxc^2*(thetac_temp(3:Nc+2)-2*thetac_temp(2:Nc+1)+thetac_temp(1:Nc))...
        - vc(it)/2/dxc*(thetac_temp(3:Nc+2)-thetac_temp(1:Nc))+ gamma.*thetan);
    
    % udpate ghost points for thetan to calculate vn
    % thetan_temp = GridUpdate_noncon(thetan,xb(it),xf(it),xb(it)-dxc,xf(it)+dxc,dxc);
    thetan_ghost_b = 2/ksigman*((eta_it+eta_st)*vb(it)*dxc*thetan(1) ...
        - eta_it*vc(it)*dxc*thetan(1) + ksigman/2*thetan(2) ...
        + (eta_it+eta_st)*dxc*Jactinb);
    thetan_ghost_f = -2/ksigman*((eta_it+eta_st)*vf(it)*dxc*thetan(Nc) ...
        - eta_it*vc(it)*dxc*thetan(Nc) - ksigman/2*thetan(Nc-1) ...
        - (eta_it+eta_st)*dxc*Jactinf);
    thetan_temp = [thetan_ghost_b; thetan; thetan_ghost_f];
    % post-process to obatin vn
    vn = eta_it/(eta_it+eta_st)*vc(it) - 1./thetan*ksigman/(eta_it+eta_st)...
        .*(thetan_temp(3:Nc+2)-thetan_temp(1:Nc))/2/dxc;
    
    %% Solute module
    Jpb = - gb*(cc(1) - cb(Nb));
    Jpf = - gf*(cc(Nc) - cf(1));
    % ghost ponints for cc
    cc_ghost_b = cc(2)    + 2*dxc/Dcc*((vb(it)-vc(it))*cc(1)  + Jpb + Jcactiveb);
    cc_ghost_f = cc(Nc-1) + 2*dxc/Dcc*((vc(it)-vf(it))*cc(Nc) + Jpf + Jcactivef);
    cc_temp = [cc_ghost_b; cc; cc_ghost_f];
    % time iteration for cc
    cc = cc + dt*(Dcc/dxc^2*(cc_temp(3:Nc+2)-2*cc_temp(2:Nc+1)+cc_temp(1:Nc)) ...
        - vc(it)/2/dxc*(cc_temp(3:Nc+2)-cc_temp(1:Nc)));
    
    % ghost ponint for cb
    cb_ghost_f = cb(Nb-1) + 2*dxb/Dcb*((vc(it)-vb(it))*cb(Nb) - Jpb - Jcactiveb);
    cb(1) = cb0;
    cb_temp = [cb; cb_ghost_f];
    % time iteraction for cb
    cb(2:Nb) = cb(2:Nb) + dt*(Dcb/dxb^2*(cb_temp(3:Nb+1)-2*cb_temp(2:Nb)+cb_temp(1:Nb-1)) ...
        - vc(it)/2/dxb*(cb_temp(3:Nb+1)-cb_temp(1:Nb-1)));
    
    % ghost point for cf
    cf_ghost_b = cf(2) + 2*dxf/Dcf*((vf(it)-vc(it))*cf(1) - Jpf - Jcactivef);
    cf(Nf) = cf0;
    cf_temp = [cf_ghost_b; cf];
    % time iteraction for cf
    cf(1:Nf-1) = cf(1:Nf-1) + dt*(Dcf/dxf^2*(cf_temp(3:Nf+1)-2*cf_temp(2:Nf)+cf_temp(1:Nf-1)) ...
        - vc(it)/2/dxf*(cf_temp(3:Nf+1)-cf_temp(1:Nf-1)));
    
    %% post-process
    ActRatio(it) = thetan(1)/thetan(Nc);
    sigma = ksigman*thetan;
    
    ActInt(it) = sum(thetan(1:Nc-1)+thetan(2:Nc))*dxc/2 ...
        + sum(thetac(1:Nc-1)+thetac(2:Nc))*dxc/2;
    ActCon(it) = ActInt(it)/(Thetat*Lc(1));
    
    %% Plot
    if mod((it-1),round(dt_plot/dt)) == 0
        iplot = iplot + 1;
        
        xcb = linspace(0,Lb(it),Nb);
        xc  = linspace(xb(it),xf(it),Nc);
        xcf = linspace(xf(it),L,Nf);
        
        figure(1) 
        subplot(4,1,1) % Solute concentration cc
        plot(xcb,cb,'k','linewidth',1); hold on
        plot(xc ,cc,'k','linewidth',1); 
        plot(xcf,cf,'k','linewidth',1);
        lim = axis;
        plot([xb(it) xb(it)],[lim(3) lim(4)],'--k','linewidth',1);
        plot([xf(it) xf(it)],[lim(3) lim(4)],'--k','linewidth',1); hold off
        title(['t =  ',num2str((it-1)*dt/60),' min'],'FontWeight','Normal');
        set(gca,'fontsize',13)
        xlabel('x (µm)','fontsize',13)
        ylabel('c (mM)','fontsize',13) 
        
        subplot(4,1,2) % Hydrostatic pressure pc
        plot([0 xb(it)],[pb0 pb_N(it)],'k','linewidth',1); hold on
        plot(xc ,pc,'k','linewidth',1); 
        plot([xf(it) L],[pf_1(it) pf0],'k','linewidth',1);
        lim = axis;
        plot([xb(it) xb(it)],[lim(3) lim(4)],'--k','linewidth',1);
        plot([xf(it) xf(it)],[lim(3) lim(4)],'--k','linewidth',1); hold off
        set(gca,'fontsize',13)
        xlabel('x (µm)','fontsize',13)
        ylabel('p_c (Pa)','fontsize',13) 
        
        subplot(4,1,3) % F-actin velocity vn
        plot(xc ,vn*1d3,'k','linewidth',1); hold off
        xlim([xb(it),xf(it)]);
        set(gca,'fontsize',13)
        xlabel('x (µm)','fontsize',13)
        ylabel('v_n (nm/s)','fontsize',13) 
        
        subplot(4,1,4) % Boundary velocities
        plot(Time/60,vb*3.6d3,'-k','linewidth',1); hold on
        plot(Time/60,vf*3.6d3,'-.k','linewidth',1); hold off
        set(gca,'fontsize',13)
        xlim([0 T_total/60])
        xlabel('time (min)','fontsize',13)
        ylabel('Velocity (µm/s)','fontsize',13) 
        legend('v_b','v_f');
        set(gcf,'color','w');
        
        % frame1(iplot) = getframe(gcf);
        
        figure(2) 
        subplot(2,1,1) % G-actin concentration thetac
        plot(xc ,thetac*1d3,'k','linewidth',1); hold off
        xlim([xb(it),xf(it)]);
        title(['t =  ',num2str((it-1)*dt/60),' min'],'FontWeight','Normal');
        set(gca,'fontsize',13)
        xlabel('x (µm)','fontsize',13)
        ylabel('\theta_c (µM)','fontsize',13) 
        
        subplot(2,1,2) % F-actin concentration thetan
        plot(xc ,thetan*1d3,'k','linewidth',1); hold off
        xlim([xb(it),xf(it)]);
        set(gca,'fontsize',13)
        xlabel('x (µm)','fontsize',13)
        ylabel('\theta_n (µM)','fontsize',13) 
        
        % frame2(iplot) = getframe(gcf);
        
        % visualize the cell with F-actin concentration coutour
        thetan_plot = [thetan'; thetan']*1d3;
        xmesh = linspace(xb(it),xf(it),Nc);
        ymesh = linspace(0,w,2);
        [Xmesh,Ymesh] = meshgrid(xmesh,ymesh);

%         figure(3)
%         subplot(2,1,1)
%         plot([0,L],[0,0],'-k','linewidth',1); hold on
%         plot([0,L],[w,w],'-k','linewidth',1); 
%         rectangle('Position',[xb(it) 0 Lc(it) w],'Curvature',w/1d2); 
%         [~,hd] = contour(Xmesh,Ymesh,thetan_plot,'fill','on'); colorbar; hold off
%         set(hd,'LevelList',linspace(min(min((thetan_plot))),max(max((thetan_plot))),25));
%         xlabel('Length (µm)')
%         ylabel('Width (µm)')
%         % set(gca,'xtick',[],'ytick',[])
%         title(['t =  ',num2str((it-1)*dt/60),' min (actin concentration in µM)'],'FontWeight','Normal');
%         set(gca,'fontsize',11)
%         box off
%         %axis equal
% 
%         frame3(iplot) = getframe(gcf);

        fprintf('Completion: %2.0f%%\n',it/Nt*1d2);

   end
end

figure(11)
subplot(2,1,1)
plot(Time/60,(vb+vf)/2*3600,'-','linewidth',1); hold off
set(gca,'fontsize',13)
xlim([0 T_total/60])
xlabel('time (min)','fontsize',13)
ylabel('Cell Velocity (µm/h)','fontsize',13)
box off
subplot(2,1,2)
plot(Time/60,ActRatio,'-','linewidth',1); hold off
set(gca,'fontsize',13)
xlim([0 T_total/60])
xlabel('time (min)','fontsize',13)
ylabel('Con. Ratio (Up/Down)','fontsize',13)
box off

figure(13)
subplot(2,1,1) % Actin Conservation factor
plot(Time/60,ActCon,'-k','linewidth',1);
set(gca,'fontsize',13)
xlim([0 T_total/60])
xlabel('time (min)','fontsize',13)
ylabel('Actin Con.','fontsize',13)

subplot(2,1,2) % Cell length ratio
plot(Time/60,Lc/Lc(1),'-k','linewidth',1);
set(gca,'fontsize',13)
xlim([0 T_total/60])
xlabel('time (min)','fontsize',13)
ylabel('Cell Length Ratio','fontsize',13)


% video = VideoWriter('moving.avi', 'Uncompressed AVI');
% video.FrameRate = 6;
% open(video)
% writeVideo(video, frame1);
% close(video);

% video = VideoWriter('actomyosin.avi', 'Uncompressed AVI');
% video.FrameRate = 6;
% open(video)
% writeVideo(video, frame2);
% close(video);

% video = VideoWriter('Color Contours.avi', 'Uncompressed AVI');
% video.FrameRate = 6;
% open(video)
% writeVideo(video, frame3);
% close(video);