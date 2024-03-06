clc,clear, close all;
%% Cerinta 1
md1 = 'penduldublu'; 

%% Cerinta 2
Tmax = 45;
t = linspace(0,Tmax,100);
u1 = 10.*double(t>=0);
usim = timeseries(u1,t);

%% Cerinta 3
load_system(md1);
set_param(md1,'StopTime',num2str(Tmax));
out = sim(md1);
simout_mat1 = out.simout_mat1;
simout_mat2 = out.simout_mat2;

simout_fcn1 = out.simout_fcn1;
simout_fcn2 = out.simout_fcn2;

% merge folosit acelasi bloc dar cu o intrare si o iesire
figure( 'Name', 'matlab function');
hold on;
plot(simout_mat1)
plot(simout_mat2)
ylabel('theta (rad)');
xlabel('Timp (s)');
legend("simout_m_a_t1" ,"simout_m_a_t2");

figure( 'Name' , 'From bloc FCN')
hold on;
plot(squeeze(simout_fcn1))
plot(squeeze(simout_fcn2))
ylabel('theta (rad)');
xlabel('Timp (s)');
legend("simout_f_c_n1" ,"simout_f_c_n2");




%% cerinta 4

err1 = norm(simout_mat1 - squeeze(simout_fcn1),2 );
err2 = norm(simout_mat2 - squeeze(simout_fcn2),2 );

%% cerinta 5 & 6


ustar = [0.1 ;0.5;0.82; 1;1.5; 1.8; 2 ; 2.8; 3.1;3.6;4;4.5;10];% luam puncte random
y1star = zeros(size(ustar));
y2star = zeros(size(ustar));

for i = 1:length(ustar)
    ustar21 = ustar(i).*double(t>=0);
    usim= timeseries(ustar21,t);
    out = sim(md1);
    y1star(i) = out.simout_mat1(end);
    y2star(i) = out.simout_mat2(end);
end    

p1 = polyfit(ustar,y1star,1); % polinomul ce trece cel mai aproape de punctele mele cele mai mici patrate
p2 = polyfit(ustar,y2star,2);
figure('Name','Aproximation from points');
plot(ustar,y1star,"*b")
hold on;
plot(ustar,y2star, "^")

ugrid = ustar(1):0.01:ustar(end);
plot(ugrid,polyval(p1,ugrid),'--r')
plot(ugrid,polyval(p2,ugrid),'--')  

xlabel("u");
ylabel("theta (radiani)");


alfa = 3;
beta = 7;
gamma = 11;
p1alfa = polyval(p1,alfa);
p2alfa = polyval(p2,alfa);
p1beta = polyval(p1,beta);
p2beta = polyval(p2,beta);
p1gamma = polyval(p1,gamma);
p2gamma = polyval(p2,gamma);



plot(alfa,polyval(p1,alfa), 'o');
hold on
plot(alfa,polyval(p2,alfa), 'o');
plot(beta,polyval(p1,beta), 'square');
plot(beta,polyval(p2,beta), 'square');
plot(gamma,polyval(p1,gamma), '*');
plot(gamma,polyval(p2,gamma), '*');
legend('(u,y1)','(u,y2)','polinom1',"polinom2",'alfa1','alfa2','beta1','beta2','gamma1','gamma2');


%% Cerinta 7
% am construit un model separat
mdl_pin = 'penduldublu_pin';

%% Cerinta 8

u0 = 17;
iu = length(u0);

[xstar,ustar,ystar,~] = trim(mdl_pin,[],u0,[],[],iu,[]);

erru = abs(ustar-u0);
%% Cerinta 9
% nu stiu ce model sa liniarizez, cel cu pini sau cel initial

[A_lin,B_lin,C_lin,D_lin] = linmod(mdl_pin, xstar, ustar);

%% Cerinta 10

vp = eig(A_lin);% este stabil deoarece are 2 perechi de valori proprii complex conjugate in semiplanul negativ
stabil = -1;
for i = 1:size(vp)
    if (real(vp(i)) > 0)
        stabil = 0;
        disp("sistem instabil");
        break;
    else
        stabil = 1;
    end
end
if(stabil == 1)
    disp("Sistem STABIL");
end

%% Cerinta 11

mdl_liniarizat = 'penduldublu_liniarizat';
r1 = 1.5.*double(t>=0);
usim = timeseries(r1,t);
load_system(mdl_liniarizat);
set_param(mdl_liniarizat,'StopTime',num2str(Tmax));
out1 = sim(mdl_liniarizat);
y_lin = out1.y_lin; %116 linii are si conditie initiala.

%% Cerinta 12


y_nl = out1.simout_mat11; %116 linii
 


err5 = norm(y_nl.Data - y_lin.Data, 'inf') 




%% cerinta 13
%u0 = 17
%[xstar,ustar,ystar,~] = trim(mdl_pin,[],u0,[],[],iu,[]);
%[A_lin,B_lin,C_lin,D_lin] = linmod(mdl_pin, xstar, ustar); liniarizarea in
%punctul u = 17

lin_model = ss(A_lin,B_lin,C_lin,D_lin); 
Te = 0.07; % discretizare tot in punctul u = 17
model_discretizat = c2d(ss(A_lin,B_lin,C_lin,D_lin),Te,'tustin'); % cu tustin stabilitatea se pastreaza dar planul stabil se muta in discul unitate pzmap(sys)
Ad = model_discretizat.A;
Bd = model_discretizat.B;
Cd = model_discretizat.C;
Dd = model_discretizat.D;
[b1,a1] = ss2tf(Ad,Bd,Cd, Dd)
% b numarator, a numitor
% linie pt a afla functia de transfer in z tf(model_discretizat)
isstable(model_discretizat)

%tf(sys_discret) adica tf(b1,a1)
% H(z) = Y(z)/U(z) =b1/a1 si de aici scot Y(z) (se da factor comun cel mai mare factor de sus si jos. -
% Y(z) = ..... 
%  1y - 3.669 z^-1Y + 5.211 z^-2Y - 3.397z^-3Y +  0.8603 z^-4Y 
% = 0.000239 * U + 2.459e-05 z^-1 U + 3.609e-15  z^-3 U + 0.0002267 z^-4 U
% 0.000239 U + 2.459e-05 z^-1 U - 0.0004411 z^-2 U  + 3.609e-15 z^-3U + 0.0002267 z^-4 U
% aplic z^-1   si o sa mai y[k] = ....  e relatia de recurenta
% pt a discretiza in simulink zero order pole
% persistent uk1 uk2 yk1 yk2 pt a pastra in memorie de-a lungul apalarii

%% Cerinta 15
mdDiscret = 'Pendul_discret'; 
load_system(mdDiscret);
set_param(mdDiscret,'StopTime',num2str(Tmax));
out = sim(mdDiscret);
simout_Liniar = out.Lin;
simout_Discret = out.Discret;
simout_Neliniar = out.Neliniar;
figure('Name','Comparison Liniar vs Neliniar vs Discret');
hold on;
plot(simout_Liniar)
plot(simout_Discret)
plot(simout_Neliniar)
ylabel('theta (rad)');
xlabel('Timp (s)');
legend("Liniar" ,"Discret","Neliniar");



