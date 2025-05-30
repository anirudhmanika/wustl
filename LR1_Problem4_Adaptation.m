clear all;close all;
vars.dt = 0.15; %utilized Euler since it was faster&more accurate for this purpose
tic

%this time is flexible in my simulation for runtime
t = [0:vars.dt:2000]; %msec
vars.t = t;
%dummy stim to populate vars.stim
inj1 = [0 0 0.5]; % amp, start, end
vars.stim = getStim(inj1,t);

vars.Cm = 1;
Vrest = -85;

%known by paper
vars.KOut = 5.4; %mM
vars.KIn = 145; %mM
vars.NaIn = 18; %mM
vars.NaOut=140; %mM
vars.CaOut=1.8; %mM
vars.CaIn=2*10^(-4); %mM
vars.PRNaK = 0.01833;

%maximal conductances
vars.gbarK = 0.282 *sqrt(vars.KOut/5.4);
vars.gbarK1 = 0.6047 *sqrt(vars.KOut/5.4);
vars.gbarNa = 23;
vars.gbarKP = 0.0183;
vars.gbarB = 0.03921;
%addition for long QT3
vars.gbarNaL = 0.01;

R = 8.3145;
F =96.485;
T = 37+273.15;
vars.const = R*T/F;

%nernst calculations
vars.ENa = (R*T/F)*log(vars.NaOut/vars.NaIn);
vars.Ek = (R*T/F)*log((vars.KOut + vars.PRNaK * vars.NaOut)/(vars.KIn + vars.PRNaK * vars.NaIn));
vars.Ek1 = (R*T/F)*log(vars.KOut/vars.KIn);
vars.Ekp = vars.Ek1;
vars.Eb = -59.87;

%define initials
[~,~,~,~,~,~, vars.minit, ~, vars.hinit, ~, vars.jinit, ~] = INAf(Vrest);
[~,~,~,~, vars.dinit, ~, vars.finit, ~] = ISIf(Vrest);
[~,~,~, vars.Xinit, ~] = IKf(Vrest);
[~, ~, vars.K1init, ~] = IK1f(Vrest,vars);

v0 = Vrest;
Ca0 = vars.CaIn;

init = [v0 Ca0];
outLQT(1,:)=init;
outCPVT(1,:)=init;

%zeta value for LQT
vars.zeta =0.5;
varsSave = vars;

numberOfHeartbeats = 10;
periods = [400:100:800];
tRecommend = [5000:1000:9000];
maxLength = max(tRecommend)./vars.dt + 1;
traces = zeros(4, length(periods), maxLength);

for j = 1:length(periods)
   period = periods(j);
   t = 0:vars.dt:tRecommend(j);
   vars.t = t;
   vars.stim = zeros(1, length(t));
   
   for i = 1:numberOfHeartbeats
       beginStim = (i-1) * period;
       inj1 = [100 beginStim beginStim+0.5];
       vars.stim = vars.stim + getStim(inj1,t);
   end
   traces(1,j,:) = [vars.stim, zeros(1, maxLength-length(vars.stim))];

   out = zeros(length(t), size(init,2));
   outLQT = out;
   outCPVT = out;
   out(1,:) = init;
   outLQT(1,:) = init;
   outCPVT(1,:) = init;
   
   % Normal simulation
   for i = 1:length(t)-1
       [dXdt, vars] = luoRudySolver(t(i), out(i,:), vars);
       out(i+1,:) = vars.dt.*dXdt + out(i,:);
   end
   vStim = out(:,1);
   traces(2,j,:) = [vStim', zeros(1, maxLength-length(vStim))];
   beginLastStim = find(t >= (numberOfHeartbeats-1)*period, 1);
   toff = t(beginLastStim:end);
   voff = vStim(beginLastStim:end);
   apdsNormal(j) = APD90(toff,voff)
   
   % LQT simulation
   for i = 1:length(t)-1
       [dXLQTdt, vars] = LR1LQT(t(i), outLQT(i,:), vars);
       outLQT(i+1,:) = vars.dt.*dXLQTdt + outLQT(i,:);
   end
   vStim = outLQT(:,1);
   traces(3,j,:) = [vStim', zeros(1, maxLength-length(vStim))];
   toff = t(beginLastStim:end);
   voff = vStim(beginLastStim:end);
   apdsLQT(j) = APD90(toff,voff)
   
   % CPVT simulation
   for i = 1:length(t)-1
       [dXLQTdt, vars] = LR1MUTATED(t(i), outCPVT(i,:), vars);
       outCPVT(i+1,:) = vars.dt.*dXLQTdt + outCPVT(i,:);
   end
   vStim = outCPVT(:,1);
   traces(4,j,:) = [vStim', zeros(1, maxLength-length(vStim))];
   toff = t(beginLastStim:end);
   voff = vStim(beginLastStim:end);
   apdsCPVT(j) = APD90(toff,voff)
end

plot(periods, apdsNormal, '-o', 'LineWidth', 1.5, 'DisplayName', 'Normal');
hold on
plot(periods, apdsLQT, '-o', 'LineWidth', 1.5,'DisplayName', 'LQTS (\zeta=0.5)');
plot(periods, apdsCPVT, '-o', 'LineWidth', 1.5, 'DisplayName', 'CPVT (0.5\cdot\beta_f)');

xlabel('Pacing Period (ms)', 'FontSize', 15);
ylabel('APD90 (ms)', 'FontSize', 15);
title('APD90 vs Pacing Period for 10 Beats', 'FontSize', 15);
legend('Location', 'northwest');
grid on;

timeElapsed = toc;

%Luo-Rudy Derivatives Output Func. (Normal Phenotype)
function [dXdt, vars] = luoRudySolver(t,x, vars)
    Iinj = interp1(vars.t, vars.stim, t, 'nearest', 'extrap');

    ddt =vars.dt;
    v=x(1);
    CaIn = x(2);
    
    [~,~,~,~,~,~, mss, taum, hss, tauh, jss, tauj] = INAf(v);
    [~,~,~,~, dss, taud, fss, tauf] = ISIf(v);
    [Xi,~,~, Xss, tauX] = IKf(v);
    [~, ~, K1ss, ~] = IK1f(v,vars);
    Kp = IKPf(v);
    
    Esi = 7.7-13.0287*log(CaIn/vars.CaOut);
    
    m = mss + (vars.minit - mss) * exp(-(ddt) / taum);
    h = hss + (vars.hinit - hss) * exp(-(ddt) / tauh);
    j = jss + (vars.jinit - jss) * exp(-(ddt) / tauj);
    d = dss + (vars.dinit - dss) * exp(-(ddt)/ taud);
    f = fss + (vars.finit - fss) * exp(-(ddt) / tauf);
    X = Xss + (vars.Xinit - Xss) * exp(-(ddt) / tauX);
    vars.minit =m;vars.hinit=h;vars.jinit=j; vars.dinit=d;vars.finit =f;vars.Xinit=X;
    
    %currents
    INa = vars.gbarNa * m^3 * h * j *(v-vars.ENa);
    Isi = 0.09 * d * f * (v-Esi);
    Ik= vars.gbarK * X * Xi * (v-vars.Ek);
    IK1 = vars.gbarK1 * K1ss * (v-vars.Ek1);
    IKP = vars.gbarKP * Kp * (v-vars.Ekp);
    Ib = vars.gbarB * (v-vars.Eb);
    
    Ii = INa + Isi + Ik + IK1 + IKP + Ib;
    
    %diff eq for voltage
    dVdt = 1/vars.Cm * (Iinj - Ii);
    
    %Calcium uptake
    dCadt = -10^(-4)*Isi+0.07*(10^(-4)- CaIn);
    
    dXdt = [dVdt, dCadt];
end

%generate stimulus function from amp, start, & end
function stim = getStim(inj,t)
    stim = zeros(size(t));
    dt = t(3)-t(2);
    start_idx = round(inj(2)/dt) + 1;
    end_idx = round(inj(3)/dt) + 1;

    stim(start_idx:end_idx) = inj(1);
end

%Ik Func.
function [Xi,aX,bX, Xss, tauX] = IKf(Vm)
    if (Vm > -100)
        Xi = 2.837 * ( exp( 0.04 * ( Vm + 77 ) ) - 1 ) / ( ( Vm + 77 ) * ( exp( 0.04 * ( Vm + 35 ) ) ) );
    else
        Xi = 1;
    end
    aX = 0.0005 * exp( 0.083 * ( Vm + 50 ) ) / ( 1 + exp( 0.057 * ( Vm + 50 ) ) );
    bX = 0.0013 * exp( -0.06 * ( Vm + 20 ) ) / ( 1 + exp( -0.04 * ( Vm + 20 ) ) );
    
    Xss = aX/(aX+bX);
    tauX = 1/(aX+bX);
end

%IK1 Func.
function [aK1, bK1, K1ss, tauK1] = IK1f(Vm, vars)
    EK1 = vars.const*log(vars.KOut/vars.KIn);

    aK1 = 1.02 / ( 1 + exp( 0.2385 * ( Vm - EK1 - 59.215 ) ) );
    bK1 = ( 0.49124 * exp ( 0.08032 * ( Vm - EK1 + 5.476 ) ) + exp ( 0.06175 * ( Vm - EK1 - 594.31 ) ) )...
        / ( 1 + exp( -0.5143 * ( Vm - EK1 + 4.753 ) ) );
    
    K1ss = aK1/(aK1+bK1);
    tauK1 = 1/(aK1+bK1);
end

%IKP func.
function [Kp] = IKPf(Vm)
    Kp = 1 / ( 1 + exp( ( 7.488 - Vm ) / 5.98) );
end

%INa Func.
function [am,bm,ah,bh,aj,bj, mss, taum, hss, tauh, jss, tauj] = INAf(Vm)
    am = 0.32*( Vm + 47.13 ) / (1 - exp( -0.1 * ( Vm + 47.13 ) ) );
    bm = 0.08 * exp( -Vm / 11 );
    if (Vm >= -40)
        ah = 0; 
        aj = 0;
        bh = 1 / ( 0.13 * ( 1 + exp( ( Vm + 10.66 ) / -11.1 ) ) );
        bj = 0.3 * exp( -2.535 * 10^-7 * Vm ) / ( 1 + exp( -0.1 * ( Vm + 32 ) ) );
    else
        ah = 0.135 * exp( ( 80 + Vm ) / -6.8 );
        bh = 3.56 * exp( 0.079 * Vm ) + 3.1 * 10^5 * exp( 0.35 * Vm );
        aj = ( -1.2714 * 10^5 * exp( 0.2444 * Vm ) - 3.474 * 10^-5 * exp( -0.04391 * Vm ) )...
            * ( Vm + 37.78 ) / ( 1 + exp( 0.311 * ( Vm + 79.23 ) ) );
        bj = 0.1212 * exp( -0.01052 * Vm ) / ( 1 + exp( -0.1378 * ( Vm + 40.14 ) ) );
    end
    
    mss = am/(am+bm);
    taum = 1/(am+bm);
    
    hss = ah/(ah+bh);
    tauh = 1/(ah+bh);
    
    jss = aj/(aj+bj);
    tauj = 1/(aj+bj);
end

%Isi func.
function [ad,bd,af,bf, dss, taud, fss, tauf] = ISIf(Vm)
    ad = 0.095 * exp( -0.01 * ( Vm - 5 ) ) / ( 1 + exp( -0.072 * ( Vm - 5 ) ) );
    bd = 0.07 * exp( -0.017 * ( Vm + 44 ) ) / ( 1 + exp( 0.05 * (Vm + 44 ) ) );
    af = 0.012 * exp( -0.008 * ( Vm + 28 ) ) / ( 1 + exp( 0.15 * ( Vm + 28 ) ) );
    bf = 0.0065 * exp( -0.02 * ( Vm + 30 ) ) / ( 1 + exp( -0.2 * ( Vm + 30 ) ) );
    
    dss = ad/(ad+bd);
    taud = 1/(ad+bd);
    
    fss = af/(af+bf);
    tauf = 1/(af+bf);
end

%PROBLEM 3, CPVT
%Luo-Rudy Derivatives Output Func. (CPVT Phenotype)
function [dXdt, vars] = LR1MUTATED(t,x, vars)
    Iinj = interp1(vars.t, vars.stim, t, 'nearest', 'extrap');

    ddt =vars.dt;
    v=x(1);
    CaIn = x(2);
    
    [~,~,~,~,~,~, mss, taum, hss, tauh, jss, tauj] = INAf(v);
    [~,~,~,~, dss, taud, fss, tauf] = ISIfMUTATED(v); %change is here
    [Xi,~,~, Xss, tauX] = IKf(v);
    [~, ~, K1ss, ~] = IK1f(v,vars);
    Kp = IKPf(v);
    
    Esi = 7.7-13.0287*log(CaIn/vars.CaOut);
    
    m = mss + (vars.minit - mss) * exp(-(ddt) / taum);
    h = hss + (vars.hinit - hss) * exp(-(ddt) / tauh);
    j = jss + (vars.jinit - jss) * exp(-(ddt) / tauj);
    d = dss + (vars.dinit - dss) * exp(-(ddt)/ taud);
    f = fss + (vars.finit - fss) * exp(-(ddt) / tauf);
    X = Xss + (vars.Xinit - Xss) * exp(-(ddt) / tauX);
    vars.minit =m;vars.hinit=h;vars.jinit=j; vars.dinit=d;vars.finit =f;vars.Xinit=X;
    
    %currents
    INa = vars.gbarNa * m^3 * h * j *(v-vars.ENa);
    Isi = 0.09 * d * f * (v-Esi);
    Ik= vars.gbarK * X * Xi * (v-vars.Ek);
    IK1 = vars.gbarK1 * K1ss * (v-vars.Ek1);
    IKP = vars.gbarKP * Kp * (v-vars.Ekp);
    Ib = vars.gbarB * (v-vars.Eb);
    INaL = vars.gbarNaL * m^3 * (v-vars.ENa);
    
    Ii = INa + Isi + Ik + IK1 + IKP + Ib + INaL;
    
    %diff eq for voltage
    dVdt = 1/vars.Cm * (Iinj - Ii);
    
    %Calcium uptake
    dCadt = -10^(-4)*Isi+0.07*(10^(-4)- CaIn);
    
    dXdt = [dVdt, dCadt];
end

%Isi Func for CPVT.
function [ad,bd,af,bf, dss, taud, fss, tauf] = ISIfMUTATED(Vm)
    ad = 0.095 * exp( -0.01 * ( Vm - 5 ) ) / ( 1 + exp( -0.072 * ( Vm - 5 ) ) );
    bd = 0.07 * exp( -0.017 * ( Vm + 44 ) ) / ( 1 + exp( 0.05 * (Vm + 44 ) ) );
    af = 1 * 0.012 * exp( -0.008 * ( Vm + 28 ) ) / ( 1 + exp( 0.15 * ( Vm + 28 ) ) );
    bf = 0.5 * 0.0065 * exp( -0.02 * ( Vm + 30 ) ) / ( 1 + exp( -0.2 * ( Vm + 30 ) ) );
    
    dss = ad/(ad+bd);
    taud = 1/(ad+bd);
    
    fss = af/(af+bf);
    tauf = 1/(af+bf);
end

%PROBLEM 2, LONG QT
%Long QT Syndrome
function [dXdt, vars] = LR1LQT(t,x, vars)
    Iinj = interp1(vars.t, vars.stim, t, 'nearest', 'extrap');

    ddt =vars.dt;
    v=x(1);
    CaIn = x(2);
    
    [~,~,~,~,~,~, mss, taum, hss, tauh, jss, tauj] = INAf(v);
    [~,~,~,~, dss, taud, fss, tauf] = ISIf(v);
    [Xi,~,~, Xss, tauX] = IKf(v);
    [~, ~, K1ss, ~] = IK1f(v,vars);
    Kp = IKPf(v);
    
    Esi = 7.7-13.0287*log(CaIn/vars.CaOut);
    
    m = mss + (vars.minit - mss) * exp(-(ddt) / taum);
    h = hss + (vars.hinit - hss) * exp(-(ddt) / tauh);
    j = jss + (vars.jinit - jss) * exp(-(ddt) / tauj);
    d = dss + (vars.dinit - dss) * exp(-(ddt)/ taud);
    f = fss + (vars.finit - fss) * exp(-(ddt) / tauf);
    X = Xss + (vars.Xinit - Xss) * exp(-(ddt) / tauX);
    vars.minit =m;vars.hinit=h;vars.jinit=j; vars.dinit=d;vars.finit =f;vars.Xinit=X;
    
    %currents
    INa = vars.gbarNa * m^3 * h * j *(v-vars.ENa);
    Isi = 0.09 * d * f * (v-Esi);
    
    %scale the delayed rectifier conductance by vars.zeta
    Ik= vars.zeta * (vars.gbarK * X * Xi) * (v-vars.Ek);
    IK1 = vars.gbarK1 * K1ss * (v-vars.Ek1);
    IKP = vars.gbarKP * Kp * (v-vars.Ekp);
    Ib = vars.gbarB * (v-vars.Eb);
    
    Ii = INa + Isi + Ik + IK1 + IKP + Ib;
    
    %diff eq for voltage
    dVdt = 1/vars.Cm * (Iinj - Ii);
    
    %Calcium uptake
    dCadt = -10^(-4)*Isi+0.07*(10^(-4)- CaIn);
    
    dXdt = [dVdt, dCadt];
end