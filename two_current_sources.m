%% Constants
eps0=8.854e-12;
mu0=4*(pi)*1e-7;
c=1/sqrt(mu0*eps0); %m/s
%%
%Create Channel-base current
I0=13e3;
eta=0.73;
tau1=-0.2e-6;
tau2=5.7e-6; %20.7e-6
n=2;
pad1=zeros(1,1000);
tf=300e-6; %bring current wave to [0,0]
N=10e3;
pad1=zeros(1,N);
Ts=tf/N;
fs=1/Ts;

t=(0:Ts:tf);
I=(I0/eta).*(((t./tau1).^n)./((t./tau1).^n+1)).*exp(-t./tau2); %Channel-base current
I=[pad1,I,pad1];

I=(I./max(I)).*12800;
N=length(I);
t=linspace(0,N*Ts,N)-tf;

plot(t*1e6,I/1000); hold all;

%% Create current at height h=300m
I0_h=I0/2;
eta_h=2*eta;
tau1_h=4*tau1;
tau2_h=2*tau2;
n=2;
pad1=zeros(1,1000);
N=10e3;
pad1=zeros(1,N);
Ts=tf/N;
fs=1/Ts;

dh=300;
v=c/3;
propagation_delay=dh/v;
prop_delay_sample=round(propagation_delay/Ts);
pad_before=[pad1 zeros(1,prop_delay_sample)];

t=(0:Ts:tf);
I_h=(I0_h/eta_h).*(((t./tau1_h).^n)./((t./tau1_h).^n+1)).*exp(-t./tau2_h); %Channel-base current
I_h=[pad_before,I_h,pad1(1:end-prop_delay_sample)];

I_h=(I_h./max(I_h)).*12800/4;
N=length(I_h);
t=linspace(0,N*Ts,N)-tf;

plot(t*1e6,I_h/1000)
xlabel('Time (us)')
ylabel('Current (kA)')
%% Constants
H=10e3; %meters
D=50e3; %meters
dz=1; %m
x0=-10;
xf=60;

if D==50e3
    dt_to_origin=467;
else if D==209e3
        dt_to_origin=996;
    else if D==250e3
            dt_to_origin=1133;
        end
    end
end

%E component before the end of the wire
pad_Dc=zeros(1,round((D/c)/Ts));
I_Dc=[pad_Dc,I]; %channel-base current retarded by D/c delay
I_h_Dc=[pad_Dc, I_h]; %current at h=300 m retarded by D/c delay

subplot(221);
time=linspace(0,length(I_Dc)*Ts,length(I_Dc)).*1e6-dt_to_origin;
plot(time,I_Dc./1e3,time,I_h_Dc./1e3);
legend('Channel-base current','Current at height h= 300 m');
xlabel('Time (\mus)');
ylabel('Current (kA)');
xlim([x0,xf]);

subplot(222);
Ez_before=-(mu0*v/(2*pi*D)).*I_Dc;
Ez_h_before=-(mu0*v/(2*pi*D)).*I_h_Dc;
Ez_before_total=Ez_before+Ez_h_before;
plot(time,-Ez_before,time,-Ez_h_before,time,-Ez_before_total);

suptitle(['D=',num2str(D/1000),' km']);

xlabel('Time (\mus)');
ylabel('Ez (V/m)');
legend('Ez before reach top of channel due to channel-base current',...
    'Ez before reach top of channel due to current at 300 m',...
    'Ez before reach top of channel due to both current sources');
xlim([x0,xf]);

Hv=round((H/v)/Ts);
pad_Hv=zeros(1,Hv);
I_Hv_Dc=[pad_Hv, I_Dc]; %channel-base current retarded by D/c+H/v delay
I_Hv_Dc=I_Hv_Dc(1:end-Hv);

I_Dc=[I_Dc, pad_Hv];
I_Dc=I_Dc(1:end-Hv);
Ez_after=-(mu0*v/(2*pi*D)).*(I_Dc-I_Hv_Dc);

I_h_Hv_Dc=[pad_Hv, I_h_Dc]; %current at h=300 m retarded by D/c+H/v delay
I_h_Hv_Dc=I_h_Hv_Dc(1:end-Hv);

I_h_Dc=[I_h_Dc, pad_Hv];
I_h_Dc=I_h_Dc(1:end-Hv);
Ez_h_after=-(mu0*v/(2*pi*D)).*(I_h_Dc-I_h_Hv_Dc);

Ez_after_total=Ez_after+Ez_h_after;

subplot(223);
plot(time,-Ez_after,time,-Ez_h_after,time,-Ez_after_total);
xlabel('Time (\mus)');
ylabel('Ez (V/m)');
legend('Ez after reach top of channel due to channel-base current',...
    'Ez after reach top of channel due to current at 300 m',...
    'Ez after reach top of channel due to both current sources');
xlim([x0,xf]);

Ez_total=Ez_before_total+Ez_after_total;
subplot(224);
plot(time,-Ez_before_total,time,-Ez_after_total,time,-Ez_total);
xlabel('Time (\mus)');
ylabel('Ez (V/m)');
legend('Ez before reach top of channel due to both current sources',...
    'Ez after reach top of channel due to both current sources',...
    'Ez total');
xlim([x0,xf]);