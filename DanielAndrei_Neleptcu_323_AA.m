%Neleptcu Daniel-Andrei 323AA Tema bonus
%༼ つ ◕_◕ ༽つ

%% Pasul 1 ⊂(◉‿◉)つ
s = tf('s');
omega = logspace(-4, 5, 1e5);

%1

P = 10/(s * (100 * s + 1) * (0.1 * s + 1));

%2

poli_P=pole(P) % Polii lui P sunt '0' '-10' '-0.01'. Deoarece nu avem poli
%in C+ sistemul este stabil. Observam ca sistemul nu este stabil in sens 
%strict pentru ca avem un pol nul in s = 0.

%3

%Avand in vedere ca avem un sistem in bucla inchisa, la alegerea uni 
%compensator static C=1 putem calcula functia de transfer a sistemului si
%functia de sensibilitate care vor conduce la constrangerea algebrica 
% S + T = 1. Mai departe, datorita prezentei unui pol cu Re(p)=0 se pot
% aplica constrangerile de interpolare si anume T(p)=1, S(p)=0.

C=1;
L=series(P,C);
S=feedback(1,L);
T=feedback(L,1);

%% Pasul 2 ⊂(◉‿◉)つ

%1

%Principiul modelului intern spune ca functia de transfer in bucla deschisa
%si anume L = P * C contine polii instabili ai referintei. Avand in vedere
%ca referinta este treapta unitara. Deoarece la formarea produsului P * C
%nu trebuie sa existe simplificari instabile in bucla de reactie, functia
%de transfer L va contine mereu polul s = '0'. In concluzie, din principiul
%modelului intern rezulta ca la alegerea oricarui regulator stabilizator se
%va realiza urmarirea asimptotica.

C=1;
L=series(P,C);
T=feedback(L,1);
figure('Name','Raspuns treapta');
step(T);

%2 
%Proiectarea lui WS
%Ne dorim eroare stationara nula la referinta treapta unitara si eroare
%<0.1% in regim stationar pentru referinte armonice in intervalul (0, 0.1)
%Conform conditiei de performanta nominala, avem |W_S * S|<1
%Prin impartirea cu |S| si constrangerea lui |S| cu eroarea de 0.1% ne
%rezulta ca |W_S|>1000


%nu suport polii din C+ ai lui P
C_1 = 250 *( (10 * s + 1) )/( (s/1000 + 1)^2);
L = series(P, C_1);
S = feedback(1, L);
T = feedback(L, 1);
W_S = 1200 /( (2 * s + 1)^4 );

a = 50;
b = 1000;
C_faza = (s/a+1)/(s/b+1); 
C_2 = series(C_1, C_faza);
L_2 = series(P, C_2);
S_2=feedback(1,L_2);
T_2=feedback(L_2,1);
W_SS=series(W_S,S_2);
conditieWS=bode(W_SS, omega);
conditieWS=squeeze(conditieWS);
conditieWS=max(conditieWS);
%ne dorim ca stabilitatea sa fie mai mica decat 1 pentru a indeplini
%conditia de performanta nominala |W_S * S|<1

figure('Name','Bode W_S');
bode(W_S, omega);
figure('Name','Nyquist L - Stabilitate bucla inchisa');
nyquist(L);

%3 

figure('Name','Bode P, W_S, C');
hold on;
bodemag(P, omega);
bodemag(W_S, omega);
hold off;
%In cazul in care avem P 'mic' si W_S va fi 'mare', vom avea functia de
%sensibilitate S 'mare', iar la conditia de performanta nominala |W_S * S|< 1
%nu va mai fi indeplinita

%% Pasul 3 ⊂(◉‿◉)つ

%1
W_T = 0.01*s/(0.001 * s+ 1);
W_TT=series(W_T,T_2);
conditieWT=bode(W_TT, omega);
conditieWT=squeeze(conditieWT);
conditieWT=max(conditieWT);
%conditia pentru |W_T * T| <1 este satisfacuta. Am realizat cand am ajuns
%la pasul 5 ca nu verificasem pentru WT si am folosit compensatorul de faza
%pentru a stabiliza aceasta conditie. Din fericire aceasta modificare a
%tinut si pentru W_S. Foarte multe modificari ale coeficientilor/numarului
%de poli le-am facut folosind in mare parte facand trial and error si dupa
%cu putin sisotool.

%2

figure('Name','Bode W_S, axa 0dB, W_T');
bodemag(W_S, omega);
hold on;
bodemag(tf(1, 1), omega);
hold on;
bodemag(W_T, omega);
%Intersectia lui W_S si W_T se afla sub axa de 0dB, deci conditia de minim
%<1 este indeplinita.

%3
%Observand diagramele bode pentru amplitudine de la subpunctul 2, pozitionand
%cursurul pe W_S si W_T se poate observa grafic ca pentru un omega mai mare
%decat 101 graficul W_T se afla deasupra axei de 0dB, iar pentru acelasi
%interval de valori, graficul W_S trece de -110dB, deci este mult mai mic
%decat axa 0dB corespunzatoare valorii 1.

%% Pasul 4 ⊂(◉‿◉)つ

omega_jf = logspace(-4, -1, 1e4);
omega_if = logspace(2, 5, 1e4);

mag_ws_jf=bode(W_S,omega_jf);
mag_ws_jf=squeeze(mag_ws_jf);

mag_wt_jf=bode(W_T,omega_jf);
mag_wt_jf=squeeze(mag_wt_jf);

R_jf = mag_ws_jf./(1-mag_wt_jf);

mag_ws_if=bode(W_S,omega_if);
mag_ws_if=squeeze(mag_ws_if);

mag_wt_if=bode(W_T,omega_if);
mag_wt_if=squeeze(mag_wt_if);

R_if = (1-mag_ws_if)./mag_wt_if;
figure('Name','Loopshaping W_S, W_T si L')
hold on;
semilogx(omega_jf,mag2db(R_jf));
semilogx(omega_if,mag2db(R_if));
bodemag(L,omega);
margin(L);
hold off;
%2
%exces poli-zerouri pt L=4 , P=3 e(L)-e(P)=1>0 good, C1 are un singur zerou
%care este negativ deci apartine in C_, din grafic se poate observa ca
%L(0)>0 si |L(jw)|>Rif si < pentru w din intervalele de joasa si inalta 
%frecventa

%3
%Transferul de la referinta la iesire si de la zgomot la iesire depind de T
%iar pentru ca T este mic la frecvente mari putem rejecta zgomotele.
%% Pasul 5 ⊂(◉‿◉)つ

%1

omega_jf = logspace(-4, -1, 1e4);
omega_if = logspace(2, 5, 1e4);

a = 50;
b = 1000;
C_faza = (s/a+1)/(s/b+1); 
C_2 = series(C_1, C_faza);
L_2 = series(P, C_2);
figure('Name','Loopshaping W_S, W_T si L2')
hold on;
semilogx(omega_jf,mag2db(R_jf));
semilogx(omega_if,mag2db(R_if));
bodemag(L_2,omega);
margin(L_2);
hold off;

S_2 = feedback(1, L_2);
T_2 = feedback(L_2, 1);
normInf = max(squeeze(bode(S_2, omega)));
margVect = 1/(normInf)

%2

% Marginea vectoriala este mai mare de decat standardul de referinta 0.5,
% utilizat in majoritatea aplicatiilor de reglare (0.79), iar pentru ca
% recomandarea este ca aceasta sa fie mai mare decat 0.5, are o valoare
% buna.In proiectare am dorit ca |S| sa se afle in jurul valorii lui 1,
% deci putem spune ca ne asteptam la o margine vectoriala in jurul valorii
% lui 1, ceea ce se poate spune ca am atins avand in vedere ca este ~0.8

%3

% Subpunctul 3 este realizat prin proiectarea compensatorului de faza,
% obtinand o margine vectoriala de ~0.8. Deci vom alege C3 in continuare
% egal cu C2
C_3=C_2;

%% Pasul 6 ⊂(◉‿◉)つ

%1
polii_S=pole(S_2)
%Toti polii lui S3=S2 au partea reala din C_, deci bucla este stabila

%2
[mag_W_SS,~]=bode(W_SS,omega);
mag_W_SS=squeeze(mag_W_SS);
[mag_W_TT,~]=bode(W_TT,omega);
mag_W_TT=squeeze(mag_W_TT);
gamma_prob=max(abs(mag_W_SS)+abs(mag_W_TT)) %(ㆆ _ ㆆ)

%3 gamma_prob < 1 si marginea vectoriala este > 0.5 deci procedura de
%obtinere a unui regulator se incheie. ʕノ•ᴥ•ʔノ ︵ ┻━┻

% ヽ༼ ຈل͜ຈ༼ ▀̿̿Ĺ̯̿̿▀̿ ̿༽Ɵ͆ل͜Ɵ͆ ༽ﾉ