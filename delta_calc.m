clear all
close all
clc

dt = 0.0001;
Tr = 0.1;
tf = 10;
t = (0:dt:(tf)).';
ti1 = 3;
ti2 = 6;
ti3 = 10;
nt = length(t);
n1 = ti1/dt;
n2 = (ti2-ti1)/dt;



a0 = 0.1;
v0 = 5;
s0 = 15;

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

noise1 = 0.0*randn(nt,1);
noise2 = 0.0*randn(nt,1);

ps_v = [];
ps_f = [];
tps = [];

for i=2:nt
    ti = dt*(i-1);
    at_s(i) = int_tustin(jt_s(i), jt_s(i-1), at_s(i-1), dt);
    vt_s(i) = int_tustin(at_s(i), at_s(i-1), vt_s(i-1), dt);
    st_s(i) = int_tustin(vt_s(i), vt_s(i-1), st_s(i-1), dt);
    % Scan Time
    if mod(ti, Tr) == 0
        [delta, ts1, ts2, grau, cx]  = pol2roots(1/2*jt_s(i), at_s(i), (vt_s(i)+ noise1(i)));
        if grau == 1
            a_ts = at_s(i) + jt_s(i)*ts1;
            if abs(a_ts) < 0.01
                ps_f = [ps_f 1];
            else
                ps_f = [ps_f 0];
            end
            ps = p23(0, 1/2*at_s(i), (vt_s(i) + noise1(i)), (st_s(i)+ noise2(i)), ts1);
            ps_v = [ps_v ps];
            tps = [tps ti];
        else
            if cx == 0
                if ts1 < 0
                    if ts2 < 0
                        ps = Inf;
                        a_ts = at_s(i) + jt_s(i)*ts2;
                    else
                        ps = p23((1/6)*jt_s(i), (1/2)*at_s(i), (vt_s(i)+noise1(i)), (st_s(i)+noise2(i)), ts2);
                        a_ts = at_s(i) + jt_s(i)*ts2;
                    end
                else
                    ps = p23((1/6)*jt_s(i), (1/2)*at_s(i), (vt_s(i)+noise1(i)), (st_s(i)+noise2(i)), ts1);
                    a_ts = at_s(i) + jt_s(i)*ts1;
                end
                if abs(a_ts) < 0.01
                    ps_f = [ps_f 1];
                else
                    ps_f = [ps_f 0];
                end
                ps_v = [ps_v ps];
                tps = [tps ti];
            end
        end
    end
end
figure
hold on
plot(t,st_s + noise2)
title('Posiçao')
xlabel('Tempo [s]')
ylabel('Posição [m]')
hold off
% 
figure
hold on
plot(t,vt_s + noise1)
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

figure
hold on
plot(tps, ps_v)
lim = axis;
plot(tps, 0.5*lim(4)*ps_f, 'r')
title('Parada')
xlabel('Tempo [s]')
ylabel('Posição [m]')
hold off
% sim('simu')
% 
