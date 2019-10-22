vm = (-1.77 .* (10.^-6))./60;
km = 28.5 .* (10.^-6);
k1 = 4.3 .* 10.^(-3);
k2 = 2.2 .* 10.^(-2); 
k3 = 2.2 .* 10.^(-3);
k4 = 1.2 .* 10.^(-2);
k5 = 4.0 .* 10.^(-4);
k6 = 3.3 .* 10.^(-1);
k7 = 4.3 .* 10.^(-3);
k8 = 9 .* 10.^(-3);
k9 = 8.3 .* 10.^(-4);
k10 = 0.05;
k11 = 2.2 .* 10.^(-3);
k12 = 2.6 .* 10.^(-2);
k13 = 4.0 .* 10.^(-4);
k14 = 0.9; 
k15 = 1.0 .* 10^(-4);
RNAP_Ecoli = 2.5 .* (10.^(-7));
promoter1 = 1.5 .* (10.^(-7));
promoter2 = 2.25 .* (10.^(-8));
tyr = 2.9 .* (10.^(-5));

syms E_out(t);
ode1 = diff(E_out, t) == (vm .* E_out)./(km + E_out);
E_outSol(t) = dsolve(ode1, E_out(0) == 5.0 .* (10 .^ (-6)));

syms E_in(t);
ode2 = diff(E_in, t) == (k1 .* E_outSol(t)) - (k1 .* E_in);
E_inSol(t) = dsolve(ode2, E_in(0) == 0);

syms mRNA_EsRNAP(t);
ode3 = diff(mRNA_EsRNAP, t) == (k2 .* RNAP_Ecoli .* promoter1) - (k3 .* mRNA_EsRNAP);
mRNA_EsRNAPSol(t) = dsolve(ode3, mRNA_EsRNAP(0) == 0);

syms RNAP(t);
ode4 = diff(RNAP, t) == (k4 .* mRNA_EsRNAPSol(t)) - (k5 .* RNAP);
RNAPSol(t) = dsolve(ode4, RNAP(0) == 0);

syms mRNA_RFP(t);
ode6 = diff(mRNA_RFP, t) == (k6 .* RNAPSol(t) .* promoter2) - (k7 .* mRNA_RFP);
mRNA_RFPSol(t) = dsolve(ode6, mRNA_RFP(0) == 0);

syms RFP(t);
ode7 = diff(RFP, t) == (k8 .* mRNA_RFPSol(t)) - (k9 .* RFP);
RFPSol(t) = dsolve(ode7, RFP(0) == 0);

syms mRNA_TAL(t);
ode10 = diff(mRNA_TAL, t) == (k10 .* RNAPSol(t) .* promoter2) - (k11 .* mRNA_TAL);
mRNA_TALSol(t) = dsolve(ode10, mRNA_TAL(0) == 0);

syms TAL(t);
ode8 = diff(TAL, t) == (k12 .* mRNA_TALSol(t)) - (k13 .* TAL);
TALSol(t) = dsolve(ode8, TAL(0) == 0);

syms pcoum(t);
ode9 = diff(pcoum, t) == (k14 .* TALSol(t) .* tyr) ./ (k15 + tyr);
pcoumSol(t) = dsolve(ode9, pcoum(0) == 0);

figure(1)
fplot(10.^(6) .* E_outSol(t), [0, 10800], 'color', 'b', 'linewidth', 1.3);
set(gca,'XTick', 0:3600:10800)
set(gca,'XTickLabel',0:1:3)
xlabel('Time (hours)', 'fontsize', 14);
ylabel('[Estrogen]_{out} (\muM)', 'fontsize', 14);

figure(2)
fplot(10.^(6) .* E_inSol(t), [0, 10800], 'color', 'r', 'linewidth', 1.3);
set(gca,'XTick', 0:3600:10800)
set(gca,'XTickLabel',0:1:3)
xlabel('Time (hours)', 'fontsize', 14);
ylabel('[Estrogen]_{in} (\muM)', 'fontsize', 14);

figure(3)
fplot(10.^(12) .* mRNA_EsRNAPSol(t), [0, 10800], 'color', 'g', 'linewidth', 1.3);
set(gca,'XTick', 0:3600:10800)
set(gca,'XTickLabel',0:1:3)
set(gca, 'YTick', 0:0.1:0.5);
xlabel('Time (hours)', 'fontsize', 14);
ylabel('[mRNA-EsRNAP] (pM)', 'fontsize', 14);

figure(4)
fplot(10.^(12) .* RNAPSol(t), [0, 21600], 'color', [1 0.5 0], 'linewidth', 1.3);
set(gca,'XTick', 0:3600:21600)
set(gca,'XTickLabel',0:1:6)
xlabel('Time (hours)', 'fontsize', 14);
ylabel('[RNAP] (pM)', 'fontsize', 14);

figure(5)
fplot(10.^(18) .* mRNA_RFPSol(t), [0, 21600], 'color', 'm', 'linewidth', 1.3);
hold on
fplot(10.^(18) .* mRNA_TALSol(t), [0, 21600], 'color', 'k', 'linewidth', 1.3);
set(gca,'XTick', 0:3600:21600)
set(gca,'XTickLabel',0:1:6)
xlabel('Time (hours)', 'fontsize', 14);
ylabel('Concentration (aM)', 'fontsize', 14);
legend('[mRNA-RFP]', '[mRNA-TAL]');

figure(6)
fplot(10.^(15) .* RFPSol(t), [0, 21600], 'color', 'c', 'linewidth', 1.3);
hold on
fplot(10.^(15) .* TALSol(t), [0, 21600], 'color', 'b', 'linewidth', 1.3);
set(gca,'XTick', 0:3600:21600)
set(gca,'XTickLabel',0:1:6)
set(gca, 'YTick', 0:0.1:0.5);
xlabel('Time (hours)', 'fontsize', 14);
ylabel('Concentration (fM)', 'fontsize', 14);
legend('[RFP]', '[TAL]');

figure(7)
fplot(10.^(12) .* pcoumSol(t), [0, 21600], 'color', 'r', 'linewidth', 1.3);
set(gca,'XTick', 0:3600:21600)
set(gca,'XTickLabel',0:1:6)
xlabel('Time (hours)', 'fontsize', 14);
ylabel('[p-coumaric acid] (pM)', 'fontsize', 14);