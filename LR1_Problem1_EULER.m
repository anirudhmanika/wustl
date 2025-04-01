clear all;close all;
%vars.dt = 0.26; %max for runge kutta
vars.dt = 0.14; %max for euler

tic

t = [0:vars.dt:600]; %msec
vars.t = t;
beginStim = 30;
inj1 = [60 beginStim beginStim+0.5]; % amp, start, end
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
[am,bm,ah,bh,aj,bj, vars.minit, taum, vars.hinit, tauh, vars.jinit, tauj] = INAf(Vrest);
[ad,bd,af,bf, vars.dinit, taud, vars.finit, tauf] = ISIf(Vrest);
[Xi,aX,bX, vars.Xinit, tauX] = IKf(Vrest);
[aK1, bK1, vars.K1init, tauK1] = IK1f(Vrest,vars);


v0 = Vrest;
Ca0 = vars.CaIn;

init = [v0 Ca0];

%Euler and midpoint methods here:
out(1,:)=init;
currents = [0 0 0 0 0 0];
iterationsLeft =length(t)-1;

for i = 1:length(t)-1
    %euler
    %need to pass through vars bc of gating updates
    [dXdt, vars, currents] = luoRudySolver(t(i), out(i,:), vars, i,currents);
        
    out(i+1,:) = vars.dt.*dXdt + out(i,:);
end
INa = currents(1,:);
Isi = currents(2,:);
Ik = currents(3,:);
IK1 = currents(4,:);
IKP = currents(5,:);
Ib = currents(6,:);

v = out(:,1);
CaIn = out(:,2);

figure;
sgtitle(['(Euler, dt=0.14) LR1 Model Response to ', num2str(max(vars.stim)), '\muA/cm^2 Injection for ',num2str(inj1(3)-inj1(2)), 'ms @ t=', num2str(inj1(2)), 'ms'], FontSize=15)

%voltage
subplot(8,1,1);
plot(t,v, LineWidth=1.5);
legend("Membrane Voltage", FontSize=13)
ylabel("Voltage (mV)")
grid on;

%ina
subplot(8,1,2);
plot(t,INa, LineWidth=1.5);
legend("I_{Na}", FontSize=13)
ylabel("Current (\muA/cm^2)")
grid on;

%isi
subplot(8,1,3);
plot(t,Isi, LineWidth=1.5);
legend("I_{si}", FontSize=13)
ylabel("Current (\muA/cm^2)")
grid on;

%ik
subplot(8,1,4);
plot(t,Ik, LineWidth=1.5);
legend("I_{K}", FontSize=13)
ylabel("Current (\muA/cm^2)")
grid on;

%ik1
subplot(8,1,5);
plot(t,IK1, LineWidth=1.5);
legend("I_{K1}", FontSize=13)
ylabel("Current (\muA/cm^2)")
grid on;

%ikp
subplot(8,1,6);
plot(t,IKP, LineWidth=1.5);
legend("I_{Kp}", FontSize=13)
ylabel("Current (\muA/cm^2)")
grid on;

%ib
subplot(8,1,7);
plot(t,Ib, LineWidth=1.5);
legend("I_b", FontSize=13)
ylabel("Current (\muA/cm^2)")
grid on;

%stimulus
subplot(8,1,8);
plot(t,vars.stim, 'Color', 'r');
ylim([0 100])
legend(['Injection of ', num2str(max(vars.stim)), '\muA/cm^2'], FontSize=13)
ylabel("Inj. Curr. (\muA/cm^2)")
grid on;

xlabel("Time (ms)", FontSize=15)

timeElapsed = toc

function stim = getStim(inj,t)
    stim = zeros(size(t));
    dt = t(3)-t(2);

    stim(inj(2)/dt + 1:inj(3)/dt +1) = inj(1);
end

%Luo-Rudy Derivatives Output Func. (Normal Phenotype)
function [dXdt, vars,currents] = luoRudySolver(t,x, vars, index,currents)
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
    INa(index) = vars.gbarNa * m^3 * h * j *(v-vars.ENa);
    Isi(index) = 0.09 * d * f * (v-Esi);
    Ik(index)= vars.gbarK * X * Xi * (v-vars.Ek);
    IK1(index) = vars.gbarK1 * K1ss * (v-vars.Ek1);
    IKP(index) = vars.gbarKP * Kp * (v-vars.Ekp);
    Ib(index) = vars.gbarB * (v-vars.Eb);
    
    Ii = INa(index) + Isi(index) + Ik(index) + IK1(index) + IKP(index) + Ib(index);
    
    %diff eq for voltage
    dVdt = 1/vars.Cm * (Iinj - Ii);
    
    %Calcium uptake
    dCadt = -10^(-4)*Isi(index)+0.07*(10^(-4)- CaIn);
    
    currents(1,index+1) = INa(index);
    currents(2,index+1) = Isi(index);
    currents(3,index+1) = Ik(index);
    currents(4,index+1) = IK1(index);
    currents(5,index+1) = IKP(index);
    currents(6,index+1) = Ib(index);
    
    dXdt = [dVdt, dCadt];
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