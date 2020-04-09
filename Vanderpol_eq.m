%% Van der Pol equation - numerical solution 
%
% y'' - mu (1 - y^2) y' + ws y = 0
%
prompt = {'Enter \mu value','Enter \omega_s value','Enter integration interval t_i:','Enter integration interval t_f:','Enter initial condition y(t_i):','Enter initial condition  y \prime (t_i):'};
dlgtitle = 'Input';
dims = [1 50; 1 50; 1 50; 1 50; 1 50; 1 50];
definput = {'2','1','0','100','2','0'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);


%temp = inputdlg('Enter \mu value');
%mu = str2double(temp{1,1});
mu = str2double(answer{1,1});
ws = str2double(answer{2,1});
%tspan = [0 100]; % integration interval [t_i t_f]
tspan = [str2double(answer{3,1}) str2double(answer{4,1})];

%y0 = [2; 0]; % initial condition [y(t_i); y'(t_i)]
y0 = [str2double(answer{5,1}); str2double(answer{6,1})];

options = odeset('MaxStep',0.5); % max integration step

ode = @(t,y) vdp_eq(t,y,mu,ws);
[t,y] = ode45(ode,tspan,y0,options);

figure(1);
plot(t,y(:,1),'b','linewidth',2); 
grid on
xlabel('$t$','Interpreter','latex','FontSize',20);
ylabel('$y$','Interpreter','latex','FontSize',20);
title(strcat('Van der Pol equation - $y(t)$ - $\mu$ = ',num2str(mu),' e $\omega_s$ = ',num2str(ws)),'Interpreter','latex','FontSize',20);

figure(2);
plot(t,y(:,2),'r','linewidth',2);
grid on
xlabel('$t$','Interpreter','latex','FontSize',20);
ylabel('$\dot y$','Interpreter','latex','FontSize',20);
title(strcat('Van der Pol equation - $\frac{dy}{dt}(t)$ - $\mu$ = ',num2str(mu),' e $\omega_s$ = ',num2str(ws)),'Interpreter','latex','FontSize',20);

figure(3);
plot(y(:,1),y(:,2));
hold on
quiver(y(:,1), y(:,2), gradient(y(:,1)), gradient(y(:,2) ));
hold off
xlabel('$y$','Interpreter','latex','FontSize',20);
ylabel('$\dot y$','Interpreter','latex','FontSize',20);
title(strcat('Van der Pol equation - $y(t)$ vs $\frac{dy}{dt}(t)$ - gradient - $\mu$ = ',num2str(mu)),'Interpreter','latex','FontSize',20);

% Fast Fourier Transform

T = 0.1;  %s % Sampling period  
Fs = 1/T; % Sampling frequency 
L = 1000; %s % Length of signal
% Time vector

tspan = 0:0.1:100; % fix step[t_i t_f]
y0 = [2; 0]; % initial condition [y(t_i); y'(t_i)]
options = odeset('MaxStep',0.5); % max integration step

ode = @(t,y) vdp_eq(t,y,mu,ws);
[t,y] = ode45(ode,tspan,y0,options);

Y = fft(y(:,1));

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure(4);
f = Fs*(0:(L/2))/L;
plot(f,P1,'linewidth',2); 
title('Single-Sided Amplitude Spectrum of $X(t)$','Interpreter','latex','FontSize',20);
xlabel('$f$ $(Hz)$','Interpreter','latex','FontSize',20);
ylabel('$| fft(y)/L |$','Interpreter','latex','FontSize',20);

kk = 1;
max = [0 0];
for ii=2:1:(length(f)-1)
    if P1(ii-1)<=P1(ii) && P1(ii+1)<=P1(ii)
        max(kk,1) = f(ii);
        max(kk,2) = P1(ii);
        kk=kk+1;
    end
end
figure(5);
for ii=1:1:length(max)
    hold on
    grid on
    plot([max(ii,1) max(ii,1)],[max(ii,2) 0],'b','linewidth',3);
end
hold off
title('Dirac Delta Function - Spectrum of $X(t)$','Interpreter','latex','FontSize',20);
xlabel('$f$ $(Hz)$','Interpreter','latex','FontSize',20);
ylabel('$| fft(y)/L |$','Interpreter','latex','FontSize',20);

function dydt = vdp_eq(~,y,mu,ws)
dydt = [y(2); mu * (1-y(1)^2)*y(2)-y(1)*ws];
end
