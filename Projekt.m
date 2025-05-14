close all
% Matematisk statistik Datorövning 2
%% Läs in data
vindelmax = readtable('vindelälven_årsmax.csv', 'VariableNamingRule', 'preserve')

luleamax = readtable('luleälven_årsmax.csv', 'VariableNamingRule', 'preserve')

%% Byt variabelnamn vindel
vindelmax.Properties.VariableNames = ["ar", "maxflode"]

%% Byt variabelnamn luleå
luleamax.Properties.VariableNames = ["ar", "maxflode"]

%% Sammanfattning vindel
summary(vindelmax)

%% Sammanfattning luleå
summary(luleamax)

%% Parameterskattningar vindel
m = mean(vindelmax.maxflode);
s = std(vindelmax.maxflode);
fprintf('Medelvärde x-streck = %.2f; standardavvikelse s =%.2f\n', m, s)
param_gumb = gumbfit(vindelmax.maxflode);
a = param_gumb(1);
b = param_gumb(2);
fprintf("Läge a = %.2f; Skala b = %.2f\n", a, b)

%% Parameterskattningar luleå
m = mean(luleamax.maxflode);
s = std(luleamax.maxflode);
fprintf('Medelvärde x-streck = %.2f; standardavvikelse s =%.2f\n', m, s)
param_gumb = gumbfit(vindelmax.maxflode);
a = param_gumb(1);
b = param_gumb(2);
fprintf("Läge a = %.2f; Skala b = %.2f\n", a, b)

%% Gumbel väntevärde och standardavvikelse
% Ersätt ? med något lämpligt
EX_gumb = 804.96 + vpa(eulergamma())*204.37;
DX_gumb = 204.37*pi/sqrt(6);
fprintf("E(X) = %.2f; D(X) = %.2f\n", EX_gumb, DX_gumb)

%% subplottar
vindelmax_p1=zeros(29, 1);
vindelmax_p2=zeros(43, 1);

for e = 1:29
vindelmax_p1(e)=vindelmax.maxflode(e);
end



for e = 1:43
vindelmax_p2(e)=vindelmax.maxflode(e+29);
end




subplot(1, 2, 1)
plot(vindelmax_p1)

subplot(1, 2, 2)
plot(vindelmax_p2)





%%

%subplottar
luleamax_p1=zeros(40, 1);
luleamax_p2=zeros(43, 1);

for e = 1:40
luleamax_p1(e)=luleamax.maxflode(e);
end



for e = 1:43
luleamax_p2(e)=luleamax.maxflode(e+40)
end




subplot(1, 2, 1)
plot(luleamax_p1)

subplot(1, 2, 2)
plot(luleamax_p2)


%%
m_vindel1 = mean(vindelmax_p1)
m_vindel2 = mean(vindelmax_p2)
m_lulea1 = mean(luleamax_p1)
m_lulea2 = mean(luleamax_p2)

s_vindel1 = std(vindelmax_p1)
s_vindel2 = std(vindelmax_p2)
s_lulea1 = std(luleamax_p1)
s_lulea2 = std(luleamax_p2)

vindel1max=max(vindelmax_p1)
vindel2max=max(vindelmax_p2)
lulea1max=max(luleamax_p1)
lulea2max=max(luleamax_p2)

vindel1min=min(vindelmax_p1)
vindel2min=min(vindelmax_p2)
lulea1min=min(luleamax_p1)
lulea2min=min(luleamax_p2)


%%

%Fråga 2:
%Hur kan man modellera  ̊arsmaximum i de b ada  ̈alvarna under de olika perioderna? 
% Kan maxflödena anses vara normalf ̈ordelade eller passar det b ̈attre med 
% Gumbelf ̈ordelning? Verkar det passa bra med samma typ av fordelning i alla 
% fyra fallen? L ̈agg ocks ̊a till de skattade parametervardena f or de anpassade 
% fordelningarna i  ̈oversiktstabellen.

param_gumb_v1 =  gumbfit(vindelmax_p1);
av1 = param_gumb_v1(1);
bv1 = param_gumb_v1(2);

param_gumb_v2=   gumbfit(vindelmax_p2);
av2 = param_gumb_v2(1);
bv2 = param_gumb_v2(2);

param_gumb_l1=   gumbfit(luleamax_p1);
al1 = param_gumb_l1(1);
bl1 = param_gumb_l1(2);

param_gumb_l2=   gumbfit(luleamax_p2);
al2 = param_gumb_l2(1);
bl2 = param_gumb_l2(2);

as=[av1 av2 al1 al2]
bs=[bv1 bv2 bl1 bl2]
ms=[m_vindel1 m_vindel2 m_lulea1 m_lulea2]
ss=[s_vindel1 s_vindel2 s_lulea1 s_lulea2]
%floden=[vindelmax_p1 vindelmax_p2 luleamax_p1 luleamax_p2]

maxfloden=[vindel1max vindel2max lulea1max lulea2max]

%for e=1:4
    x = linspace(0, maxfloden(1));
% weibull
p_norms = normpdf(x, ms(1), ss(1));
F_norms = normcdf(x, ms(1), ss(1));
% Gumbel
p_gumbs = gumbpdf(x, as(1), bs(1));
F_gumbs = gumbcdf(x, as(1), bs(1));
subplot(2, 2, 1)
hold on
histogram(vindelmax_p1, 'Normalization', 'pdf')
plot(x, p_norms, 'r', 'Linewidth', 1)
plot(x, p_gumbs, 'b', 'Linewidth', 1)
legend('norm', 'gumb');
hold off

subplot(2, 2, 2)
hold on
cdfplot(vindelmax_p1)
plot(x, F_norms, 'r', 'Linewidth', 1)
plot(x, F_gumbs, 'b', 'Linewidth', 1)
hold off
subplot(2, 2, 3)
probplot('normal', vindelmax_p1)
grid on
subplot(2, 2, 4)
gumbplot(vindelmax_p1)
grid on

%end
%%


   x = linspace(0, maxfloden(2));
% weibull
p_norms = normpdf(x, ms(2), ss(2));
F_norms = normcdf(x, ms(2), ss(2));
% Gumbel
p_gumbs = gumbpdf(x, as(2), bs(2));
F_gumbs = gumbcdf(x, as(2), bs(2));
subplot(2, 2, 1)
hold on
histogram(vindelmax_p2, 'Normalization', 'pdf')
plot(x, p_norms, 'r', 'Linewidth', 1)
plot(x, p_gumbs, 'b', 'Linewidth', 1)
legend('norm', 'gumb');
hold off

subplot(2, 2, 2)
hold on
cdfplot(vindelmax_p2)
plot(x, F_norms, 'r', 'Linewidth', 1)
plot(x, F_gumbs, 'b', 'Linewidth', 1)
hold off
subplot(2, 2, 3)
probplot('normal', vindelmax_p2)
grid on
subplot(2, 2, 4)
gumbplot(vindelmax_p2)
grid on

%%
 x = linspace(0, maxfloden(3));
% weibull
p_norms = normpdf(x, ms(3), ss(3));
F_norms = normcdf(x, ms(3), ss(3));
% Gumbel
p_gumbs = gumbpdf(x, as(3), bs(3));
F_gumbs = gumbcdf(x, as(3), bs(3));
subplot(2, 2, 1)
hold on
histogram(luleamax_p1, 'Normalization', 'pdf')
plot(x, p_norms, 'r', 'Linewidth', 1)
plot(x, p_gumbs, 'b', 'Linewidth', 1)
legend('Data', 'Norm', 'Gumb');
hold off

subplot(2, 2, 2)
hold on
cdfplot(luleamax_p1)
plot(x, F_norms, 'r', 'Linewidth', 1)
plot(x, F_gumbs, 'b', 'Linewidth', 1)
hold off
subplot(2, 2, 3)
probplot('normal', luleamax_p1)
grid on
subplot(2, 2, 4)
gumbplot(luleamax_p1)
grid on

%%
 x = linspace(0, maxfloden(4));
% weibull
p_norms = normpdf(x, ms(4), ss(4));
F_norms = normcdf(x, ms(4), ss(4));
% Gumbel
p_gumbs = gumbpdf(x, as(4), bs(4));
F_gumbs = gumbcdf(x, as(4), bs(4));
subplot(2, 2, 1)
hold on
histogram(luleamax_p2, 'Normalization', 'pdf')
plot(x, p_norms, 'r', 'Linewidth', 1)
plot(x, p_gumbs, 'b', 'Linewidth', 1)
legend('norm', 'gumb');
hold off

subplot(2, 2, 2)
hold on
cdfplot(luleamax_p2)
plot(x, F_norms, 'r', 'Linewidth', 1)
plot(x, F_gumbs, 'b', 'Linewidth', 1)
hold off
subplot(2, 2, 3)
probplot('normal', luleamax_p2)
grid on
subplot(2, 2, 4)
gumbplot(luleamax_p2)
grid on


%%
% %Fråga 3. Har förvantade  ̊arliga maximala flödet i Vindel ̈alven respektive i Luleälven  ̈andrat sig mellan
%perioderna? Vi vill ha både skattning av förändringen, konfidensintervall och test. Motivera
%också om det verkar rimligt att anta att standardavvikelserna kan vara samma under båda period-
%erna f̈ör respektive  ̈alv och varför det går att använda normalapproximation.
medelskillnad_lulea = m_lulea2-m_lulea1
medelskillnad_vindel = m_vindel2 - m_vindel1

n1 = length(luleamax_p1);
n2 = length(luleamax_p2);
n3 = length(vindelmax_p1);
n4 = length(vindelmax_p2);

s_p_vindel = sqrt(( (n3-1) * s_vindel1^2 + (n4-1) * s_vindel2^2 ) / (n3 + n4 -2))
d_vindel = s_p_vindel * sqrt(1/n3 + 1/n4);
f_vindel = n3 + n4-2

f_lule = ((s_lulea1^2 / n1 + s_lulea2^2 / n2)^2) / ( (s_lulea1^2/n1)^2 / (n1-1)   +  (s_lulea2^2/n2)^2 / (n2-1)   )
d_lule = sqrt(s_lulea1^2/n1 + s_lulea2^2/n2)

alfa = 0.05/2;
tkvantil_vindel = tinv(1-alfa, f_vindel)
tkvantil_lule = tinv(1-alfa, f_lule)

I_lagre_vindel = medelskillnad_vindel - tkvantil_vindel * d_vindel
I_hogre_vindel = medelskillnad_vindel + tkvantil_vindel * d_vindel

I_lagre_lule = medelskillnad_lulea - tkvantil_lule * d_lule
I_hogre_lule = medelskillnad_lulea + tkvantil_lule * d_lule


[h, p, ci, stats] = ttest2(luleamax_p2, luleamax_p1, "Vartype","unequal");
[h, p, ci, stats] = ttest2(vindelmax_p2, vindelmax_p1, "Vartype","equal");

%%
%4. Lule ̈alvens avrinningsomr ̊ade  ̈ar dubbelt s ̊a stort som Vindel ̈alvens s ̊a man kan f ̈orv ̈anta sig att
%h ̈oga fl ̈oden i Lule ̈alven  ̈ar ungef ̈ar dubbelt s ̊a h ̈oga som h ̈oga fl ̈oden i Vindel ̈alven. Vi definierar
%”h ̈ogt fl ̈ode” som ett  ̊arsmaximum  ̈over 2000 m3/s i Lule ̈alven respektive  ̈over 1000 m3/s i Vin-
%del ̈alven. Under perioden f ̈ore regleringen (–1939) var sannolikheten f ̈or ”h ̈ogt fl ̈ode” ungef ̈ar
%30 % i b ̊ada  ̈alvarna.
%Under hur m ̊anga av  ̊aren efter reglering (1980–) var det ”h ̈ogt fl ̈ode” i Lule ̈alven respektive
%i Vindel ̈alven? Har det skett en signifikant minskning av sannolikheten f ̈or ”h ̈ogt fl ̈ode” i
%Lule ̈alven? I Vindel ̈alven?

lulea_efter = sum(luleamax_p2 >= 2000)
andel_lulea_efter = lulea_efter / p2

vindel_efter = sum(vindelmax_p2 > 1000)
andel_vindel_efter = vindel_efter / n2

binocdf(12,43,0.3)
binocdf(1,43,0.3)
