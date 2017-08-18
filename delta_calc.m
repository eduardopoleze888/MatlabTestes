clear all
close all
clc

eps = 0.2;       % Período de atualização (período de amostragem, task rate) [s];

dt = 0.001;
Tr = 0.1;
tf = 30;
t = (0:dt:(tf)).';
nt = length(t);
n1 = round(nt/3);
n2 = n1;

ti1 = 10;
ti2 = 20;
ti3 = 30;

a0 = 4;
v0 = 10;
s0 = 0;

al = -(2*v0 + a0*ti1)/(ti2+ti3-ti1);


j1 = (al-a0)/ti1;
j2 = 0;
j3 = -al/(ti3-ti2);


jt_s = [j1*ones(n1,1); j2*ones(n2,1); j3*ones(nt-n1-n2,1)];
jt_in.signals.values = [j1*ones(n1,1); j2*ones(n2,1); j3*ones(nt-n1-n2,1)];
jt_in.time = t;



st_s = zeros(nt,1);
vt_s = zeros(nt,1);
at_s = zeros(nt,1);

st_s(1) = s0;
vt_s(1) = v0;
at_s(1) = a0;

ps1_v = [];
ps2_v = [];
tps = [];

for i=2:nt
    ti = dt*(i-1);
    at_s(i) = int_tustin(jt_s(i), jt_s(i-1), at_s(i-1), dt);
    vt_s(i) = int_tustin(at_s(i), at_s(i-1), vt_s(i-1), dt);
    st_s(i) = int_tustin(vt_s(i), vt_s(i-1), st_s(i-1), dt);
    if mod(ti, Tr) == 0
        if jt_s(i) == 0
            ts1 = -vt_s(i)/at_s(i);
            ts2 = ts1;
            ps1 = st_s(i) + vt_s(i)*ts1 + 1/2*at_s(i)*ts1^2;
            ps2 = ps1;
            ps1_v = [ps1_v ps1];
            ps2_v = [ps2_v ps2];
            tps = [tps ti];
        else
            delta = at_s(i)^2 - 2*vt_s(i)*jt_s(i);
            if delta >= 0
                ts1 = (-at_s(i) + sqrt(delta))/jt_s(i);
                ps1 = st_s(i) + vt_s(i)*ts1 + (1/2)*at_s(i)*ts1^2 + (1/6)*jt_s(i)*ts1^3
                ps1_v = [ps1_v ps1];
                ts2 = (-at_s(i) - sqrt(delta))/jt_s(i);
                ps2 = st_s(i) + vt_s(i)*ts2 + (1/2)*at_s(i)*ts2^2 + (1/6)*jt_s(i)*ts2^3
                ps2_v = [ps2_v ps2];
                tps = [tps ti];
            end
        end
    end  
end

figure
hold on
plot(tps, ps1_v, tps, ps2_v)
title('Posiçao')
legend('ps1', 'ps2')
xlabel('Tempo [s]')
ylabel('Posição [m]')
hold off
% sim('simu')
% 
figure
hold on
plot(t,st_s)
title('Posiçao')
xlabel('Tempo [s]')
ylabel('Posição [m]')
hold off
% 
figure
hold on
plot(t,vt_s)
title('Velocidade')
xlabel('Tempo [s]')
ylabel('Velocidade [m/s]')
hold off
% 
figure
hold on
plot(t,at_s)
title('Aceleração')
xlabel('Tempo [s]')
ylabel('Aceleração [m/s^2]')
hold off

figure
hold on
plot(t,jt_s)
title('Jerk')
xlabel('Tempo [s]')
ylabel('Jerk [m/s^3]')
hold off
