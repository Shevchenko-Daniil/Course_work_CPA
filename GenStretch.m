clear; clc;

Tin = 50*10^(-15);                                  %длительность сгенерированного импульса
Tchirp = 100*10^(-12);                              %длительность чирпированного импульса
lambda = 800*10^(-9);                                %длина волны
c = 3*10^(8);                                       %скорость света
H = 1200000;                                        %число штрихов на метр

a = ((Tin)^2)/(8*log(2));                           %характеризует начальную длительность
b = -sqrt( a*(Tchirp)^2/(8*log(2)) - a^2);           %задает чирп для длительности Tchirp пс
w0 = 2*pi*c/lambda;                                 %центральная угловая частота
N = 2^17;                                           %количесвто точек сетки

gamma = 43*pi/180; %9.969                                  %угол падения

G = 0.3065;             %0.65                                          %расстояние между решетками

% Спектр фемтосекундного импульса

%w = linspace(0.85*w0, 1.2*w0, N);

w = 0.0*w0:(0.09*w0/N):8.0*w0;
N = length(w);
X = exp(-0.5*a.*(w).^2);

index_w_min = find( abs(w-2.1e15) <= (max(w) - min(w))/(2*N) );
index_w_max = find( abs(w-2.6e15) <= (max(w) - min(w))/(2*N) );

wleft = 2.1e15;
wright = 2.6e15;

figure(1)
    subplot(2, 1, 1);
    plot(w, X)
    title('Cпектр фемтосекундного импульса')
    xlim([wleft wright])
    subplot(2, 1, 2);
    plot(w, unwrap(angle(X)))
    title('Фаза фемтосекундного импульса')
    xlim([wleft wright])
    
    
%% Фемтосекундный импульс

dt = 2*pi/(max(w)-min(w));
t = ((1:N) - N/2 - 1)*dt;

S = real(fftshift(fft(X)));

tleft = -8.0*10^(-14);
tright = 8.0*10^(-14);

figure(2)
    plot(t, S/max(S))
    title("Начальный импульс")
    xlim([tleft tright])
    
Sin = (abs(S)).^2;
index_Sin_max = find(Sin == max(Sin));
Sin = Sin((index_Sin_max - 30000):(index_Sin_max + 30000));
    
    
    
 %% Стретчер
 
sinT = H*lambda - sin(gamma);
cosT = ( 1 - (sinT).^2 ).^0.5;
%tgT = sinT./cosT;



phi_gen_stretch = (2*G/c)*(cos(gamma).*w + ((w.^2).*cos(gamma).^2 - (2*pi*c*H).^2 + 2.*2.*pi.*c.*H.*(sin(gamma)).*(w)).^0.5 );

tmp_phi_str = zeros(size(phi_gen_stretch));
tmp_phi_str(index_w_min:index_w_max) = phi_gen_stretch(index_w_min:index_w_max);
phi_gen_stretch = tmp_phi_str;

Xch = exp(-0.5*a.*(w-w0).^2 + 1i*phi_gen_stretch);
%Xch = X.*exp(-0.5*i*b.*(w-w0).^2);

figure(3)
    subplot(2, 1, 1);
    plot(w, real(Xch))
    title('Cпектр растянутого импульса')
    xlim([wleft wright])
    subplot(2, 1, 2);
    plot(w, unwrap(angle(Xch)))
    title('Фаза растянутого импульса')
    xlim([wleft wright])
    
    
    
S = real(fftshift(ifft(Xch)));

dt1 = 2*pi/(max(w) - min(w));
t1 = ((1:N) - N/2 - 1)*dt1;
    
figure(4)
    plot(t1, (S/max(S)))
    title("Растянутый импульс через общ. формулу")
    %xlim([tleft tright])
    
    
%%
%Glist = 1:1:20;
%%

gamma = 44*pi/180;

G0 = 0.3065; %0.3065->0.3095->0.3099
for k = 1:20

G=G0 - 0.01 + k*0.001;
Glist(k) = G;

sinT = H*lambda - sin(gamma);
cosT = ( 1 - (sinT).^2 ).^0.5;
tgT = sinT./cosT;

phi_gen = (2*G/c)*(cos(gamma).*w + ((w.^2).*cos(gamma).^2 - (2*pi*c*H).^2 + 2.*2.*pi.*c.*H.*(sin(gamma)).*(w)).^0.5 );

tmp_phi = zeros(size(phi_gen));
tmp_phi(index_w_min:index_w_max) = phi_gen(index_w_min:index_w_max);
phi_gen = tmp_phi;

%Xcom = Xch.*exp(1i*phi_gen);
Xcom = exp(-0.5*a.*(w-w0).^2 + 1i*(phi_gen_stretch-phi_gen));

% figure(3)
%     subplot(2, 1, 1);
%     plot(w, real(Xcom))
%     title('Cпектр сжатого импульса через общ. формулу')
%     xlim([wleft wright])
%     subplot(2, 1, 2);
%     plot(w, unwrap(angle(Xcom)))
%     title('Фаза сжатого импульса через общ. формулу')
%     xlim([wleft wright])

    
S = (fftshift(ifft(Xcom)));



dt1 = 2*pi/(max(w) - min(w));
t1 = ((1:N) - N/2 - 1)*dt1;
    
% figure(4)
%     plot(t1, (S/max(S)))
%     title("Сжатый импульс через общ. формулу")

    
    
smax = max(abs(S));

halfMax = smax/2;
index1 = find(abs(S) >= halfMax, 1, 'first');
index2 = find(abs(S) >= halfMax, 1, 'last');

FWHMout_gen(k) = t1(index2) - t1(index1);
end

%%
% figure(50)
% plot(Glist, FWHMout_gen)

%%


index_T_min = find(FWHMout_gen == min(FWHMout_gen));
G0 = Glist(index_T_min);
%

for k = 1:20

G=G0 - 0.001 + k*0.0001;
Glist(k) = G;

sinT = H*lambda - sin(gamma);
cosT = ( 1 - (sinT).^2 ).^0.5;
tgT = sinT./cosT;

phi_gen = (2*G/c)*(cos(gamma).*w + ((w.^2).*cos(gamma).^2 - (2*pi*c*H).^2 + 2.*2.*pi.*c.*H.*(sin(gamma)).*(w)).^0.5 );

tmp_phi = zeros(size(phi_gen));
tmp_phi(index_w_min:index_w_max) = phi_gen(index_w_min:index_w_max);
phi_gen = tmp_phi;

%Xcom = Xch.*exp(1i*phi_gen);
Xcom = exp(-0.5*a.*(w-w0).^2 + 1i*(phi_gen_stretch-phi_gen));

% figure(3)
%     subplot(2, 1, 1);
%     plot(w, real(Xcom))
%     title('Cпектр сжатого импульса через общ. формулу')
%     xlim([wleft wright])
%     subplot(2, 1, 2);
%     plot(w, unwrap(angle(Xcom)))
%     title('Фаза сжатого импульса через общ. формулу')
%     xlim([wleft wright])

    
S = (fftshift(ifft(Xcom)));



dt1 = 2*pi/(max(w) - min(w));
t1 = ((1:N) - N/2 - 1)*dt1;
    
% figure(4)
%     plot(t1, (S/max(S)))
%     title("Сжатый импульс через общ. формулу")

    
    
smax = max(abs(S));

halfMax = smax/2;
index1 = find(abs(S) >= halfMax, 1, 'first');
index2 = find(abs(S) >= halfMax, 1, 'last');

FWHMout_gen(k) = t1(index2) - t1(index1);
end
%%
figure(50)
plot(Glist, FWHMout_gen, 'LineWidth', 2)
ylabel('T_{out}, c', 'FontSize', 13)
xlabel('G, м', 'FontSize', 13)
grid on

%%

index_T_min = find(FWHMout_gen == min(FWHMout_gen));
G0 = max(Glist(index_T_min));

%

for k = 1:20

G=G0 - 0.0001 + k*0.00001;
Glist(k) = G;

sinT = H*lambda - sin(gamma);
cosT = ( 1 - (sinT).^2 ).^0.5;
tgT = sinT./cosT;

phi_gen = (2*G/c)*(cos(gamma).*w + ((w.^2).*cos(gamma).^2 - (2*pi*c*H).^2 + 2.*2.*pi.*c.*H.*(sin(gamma)).*(w)).^0.5 );

tmp_phi = zeros(size(phi_gen));
tmp_phi(index_w_min:index_w_max) = phi_gen(index_w_min:index_w_max);
phi_gen = tmp_phi;

%Xcom = Xch.*exp(1i*phi_gen);
Xcom = exp(-0.5*a.*(w-w0).^2 + 1i*(phi_gen_stretch-phi_gen));

% figure(3)
%     subplot(2, 1, 1);
%     plot(w, real(Xcom))
%     title('Cпектр сжатого импульса через общ. формулу')
%     xlim([wleft wright])
%     subplot(2, 1, 2);
%     plot(w, unwrap(angle(Xcom)))
%     title('Фаза сжатого импульса через общ. формулу')
%     xlim([wleft wright])

    
S = (fftshift(ifft(Xcom)));



dt1 = 2*pi/(max(w) - min(w));
t1 = ((1:N) - N/2 - 1)*dt1;
    
% figure(4)
%     plot(t1, (S/max(S)))
%     title("Сжатый импульс через общ. формулу")

    
    
smax = max(abs(S));

halfMax = smax/2;
index1 = find(abs(S) >= halfMax, 1, 'first');
index2 = find(abs(S) >= halfMax, 1, 'last');

FWHMout_gen(k) = t1(index2) - t1(index1);
end


index_T_min = find(FWHMout_gen == min(FWHMout_gen));
G = max(Glist(index_T_min));

FWHMout = min(FWHMout_gen);


%

% G1=(G0 - 0.0001 + 0.00001):0.00001:(G0 + 0.0001);
% figure(10)
% plot(G1, FWHMout_gen)

%
%G = 0.30993; 

sinT = H*lambda - sin(gamma);
cosT = ( 1 - (sinT).^2 ).^0.5;
tgT = sinT./cosT;

phi_gen = (2*G/c)*(cos(gamma).*w + ((w.^2).*cos(gamma).^2 - (2*pi*c*H).^2 + 2.*2.*pi.*c.*H.*(sin(gamma)).*(w)).^0.5 );

tmp_phi = zeros(size(phi_gen));
tmp_phi(index_w_min:index_w_max) = phi_gen(index_w_min:index_w_max);
phi_gen = tmp_phi;

%Xcom = Xch.*exp(1i*phi_gen);
Xcom = exp(-0.5*a.*(w-w0).^2 + 1i*(phi_gen_stretch-phi_gen));

% figure(3)
%     subplot(2, 1, 1);
%     plot(w, real(Xcom))
%     title('Cпектр сжатого импульса через общ. формулу')
%     xlim([wleft wright])
%     subplot(2, 1, 2);
%     plot(w, unwrap(angle(Xcom)))
%     title('Фаза сжатого импульса через общ. формулу')
%     xlim([wleft wright])

    
S = (fftshift(ifft(Xcom)));

%pulses(k) = S;

dt1 = 2*pi/(max(w) - min(w));
t1 = ((1:N) - N/2 - 1)*dt1;

index_S_max = find((abs(S)).^2 == max((abs(S)).^2));

S8 = (abs(S)).^2;
FWHMout8 = FWHMout;

% figure(4)
%     plot(t1, (abs(S/max(S))).^2)
%     title("Сжатый импульс через общ. формулу")
%     %xlim([-12e-14 12e-14])

%%
%index_S1_max = find(S1 == max(S1));
%index_S2_max = find(S2 == max(S2));
%index_S3_max = find(S3 == max(S3));
%index_S4_max = find(S4 == max(S4));
index_S5_max = find(S5 == max(S5));
index_S6_max = find(S6 == max(S6));
index_S7_max = find(S7 == max(S7));
index_S8_max = find(S8 == max(S8));

%S1 = S1((index_S1_max - 30000):(index_S1_max + 30000));
%S2 = S2((index_S2_max - 30000):(index_S2_max + 30000));
%S3 = S3((index_S3_max - 30000):(index_S3_max + 30000));
%S4 = S4((index_S4_max - 30000):(index_S4_max + 30000));
S5 = S5((index_S5_max - 30000):(index_S5_max + 30000));
S6 = S6((index_S6_max - 30000):(index_S6_max + 30000));
S7 = S7((index_S7_max - 30000):(index_S7_max + 30000));
S8 = S8((index_S8_max - 30000):(index_S8_max + 30000));

%%
FWHMlist = [FWHMout8, FWHMout6, FWHMout4, FWHMout3, Tin, FWHMout2, FWHMout6, FWHMout5, FWHMout7];
anglelist = 41:0.5:45;

figure(20)
    plot(anglelist, FWHMlist)
    ylabel('T_{out}, фс', 'FontSize', 13)
    xlabel('\gamma°', 'FontSize', 13)
    grid on
    
    set(findall(figure(2),'type','axes'),'fontsize',16)
    %%
figure(80)   
yyaxis left
plot(S1/max(S1))
hold on
yyaxis right
plot(S5/max(S5))
hold on
plot(Sin/max(Sin))
hold off






 