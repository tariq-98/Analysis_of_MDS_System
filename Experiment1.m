&& Author@Tariq Al-Radhi

clc
clear 
close all
%% Find the Spring Constant K (Newtons per meter)
% Initialisation
load = [0 5 10,15 20 25];
voltage = [0 92.44 184.1 281.3 367.4 464.7];
volToMet = 1/350;
AB=125;
AC=690;
% Processing
BEmm = voltage*volToMet;
CDmm = (BEmm*AC)/AB;
CDm = CDmm*1e-3;
k = (load(1)-load(6))/(CDm(1)-CDm(6));
% Post Processing
figure(1)
plot(CDm,load,'-Ok','LineWidth',2);
title('Graph of Load vs Extension');
xlabel('Extension(m)');
ylabel('Load(N)');
LEtable = table(load',voltage',CDmm',CDm','VariableNames',{'Load (N)',...
    'Volatage (mv)','Extension (mm)','Extension (m)'});
disp(LEtable)
%% Time Period to find Natural Frequency (rad/sec)
% Initialisation
tp = [139.6 140.2 137.2 135.3 137.2];
tpSum = 0;
% Processing
for i = 1:length(tp)
    tpSum = tpSum + tp(i);
end
tpsAvg = tpSum/length(tp);
nf = (1/tpsAvg)*1e3*2*pi;
%% Mass of System (Kg)
mass = k/nf^2;
%% Critical Damping Coefficient
bcr = 2*sqrt(mass*k);
%% Frequency against Amplitude without Damper
% Iniitialisation
forcedFrequency = [2.51 3.01 3.5 4 4.5 5 5.5 6 6.5 6.7 6.9 ...
    7.1 7.3 7.425 7.5 8 8.5 9 9.5];
Amplitude = [300 360.6 380.8 350 400 430 430 442 526 644 ...
    852 2232 6500 6800 6876 694 530 509.8 443.1];
angularVelocity = forcedFrequency*2*pi;
% Processing
halfpower = max(Amplitude)/sqrt(2);
%solving for w1(omega 1)
x1 = zeros(1,6);
y1 = zeros(1,6);
for i = 10:15
    x1(i-9) = angularVelocity(i);
    y1(i-9) = Amplitude(i);
end
wmax = interp1(y1,x1,max(Amplitude));
w1 = interp1(y1,x1,halfpower);
%solving for w2(omega 2)
x2 = zeros(1,4);
y2 = zeros(1,4);
for i =15:19
    x2(i-14) = angularVelocity(i);
    y2(i-14) = Amplitude(i);
end
w2 = interp1(y2,x2,halfpower);
dRatio = (w2-w1)/(2*wmax);
% Damping Coefficient
b = dRatio * bcr;
% Post Processing
%ploting
figure(2)
plot(angularVelocity, Amplitude,'-ok','LineWidth',2)
hold on
line([w1,w2],[halfpower,halfpower],'LineStyle','--','LineWidth',1,...
    'color','k')
hold on 
line([w1,w1],[0,halfpower],'LineStyle','--','LineWidth',1,...
    'color','k')
hold on
line([w2,w2],[0,halfpower],'LineStyle','--','LineWidth',1,...
    'color','k')
hold on
plot(w1,halfpower,'xk')
hold on
plot(w2,halfpower,'xk')
hold on
plot(wmax,max(Amplitude),'xk')
xlabel('Angular Velocity (rad/s)')
ylabel('Amplitude (mV)')
title('Amplitude vs Angular Velocity without Damper')
legend('Amplitude vs Angular Velocity','Bandwidth')

UndampedDataTable = table(wmax,mass,k,b,bcr,dRatio,'VariableNames'...
    ,{'Natural Frequency (rad/sec)','Mass (Kg)','Spring Stiffness (N/m)'...
    ,'Damping Coefficient (Ns/m)','Critial Damping coefficition (Ns/m)'...
    ,'Damping Ratio'});
disp("Results for Undamped System")
disp(UndampedDataTable)

% Frequency against Amplitude with Damper
%% Initialisation
forcedFrequency = [2.5 3 3.5 4 4.5 5 5.5 6 6.5 6.7 6.8 7.06 ...
    7.2 7.5 7.7 8 8.5 9.3];
Amplitude = [345.7 368.5 412.4 345.7 423.8 386.1 411.5 479.1 ...
    859.1 1608 3257 1313 990.6 676.5 551.9 531.7 435.2 456.3];
angularVelocity = forcedFrequency*2*pi;
phaseDelay = Amplitude*0;
%% Processing
HalfPower = max(Amplitude)/sqrt(2);
%solving for w1(omega 1)
x1 = zeros(1,6);
y1 = zeros(1,6);
for i = 6:11
    x1(i-5) = angularVelocity(i);
    y1(i-5) = Amplitude(i);
end
wdmax = interp1(y1,x1,max(Amplitude));
w1 = interp1(y1,x1,HalfPower);
%solving for w2(omega 2)
x2 = zeros(1,3);
y2 = zeros(1,3);
for i =11:13
    x2(i-10) = angularVelocity(i);
    y2(i-10) = Amplitude(i);
end
w2 = interp1(y2,x2,HalfPower);
dRatio2 = (w2-w1)/(2*wdmax);
% Finding Phase Delay
for i=1:length(phaseDelay)
    phaseDelay(i) = atan(-1*(2*dRatio2*(angularVelocity(i)/wdmax)) ... 
        /(1-(angularVelocity(i)/wdmax)^2));
     if(phaseDelay(i)>0)
        phaseDelay(i) = phaseDelay(i) - pi;
     end
end


%% Post Processing
%plotting angular velocity against amplitude
figure(3)
plot(angularVelocity,Amplitude,'-ok','LineWidth',2)
hold on
line([w1,w2],[HalfPower,HalfPower],'LineStyle','--','LineWidth',1,...
    'color','k')
hold on 
line([w1,w1],[0,HalfPower],'LineStyle','--','LineWidth',1,...
    'color','k')
hold on
line([w2,w2],[0,HalfPower],'LineStyle','--','LineWidth',1,...
    'color','k')
hold on
plot(w2,HalfPower,'xk')
hold on
plot(w1,HalfPower,'xk')
hold on
plot(wdmax,max(Amplitude),'xr')
hold on
xlabel('Angular Velocity (rad/s)')
ylabel('Amplitude (mV)')
title('Amplitude vs Angular Velocity with Damper')
legend('Amplitude (mv) Against Angular Velocity (rad/sec)','Bandwidth line')


% Plotting phase against angular velocity
figure(4)
plot(angularVelocity,phaseDelay,'-ok','LineWidth',2)
xlabel('Angular Velocity (rad/s)')
ylabel('Phase Delay (rad)')
title('Phase Delay vs Angular Velocity with Damper')
mass = k/wdmax^2;
bcr = 2*sqrt(mass*k);
b = dRatio2*bcr;
disp("Results for Damped System")
DampedDataTable = table(wdmax,mass,k,b,bcr,dRatio2,'VariableNames'...
    ,{'Natural Frequency (rad/sec)','Mass (Kg)','Spring Stiffness (N/m)'...
    ,'Damping Coefficient (Ns/m)','Critial Damping coefficition (Ns/m)'...
    ,'Damping Ratio'});
disp(DampedDataTable)

