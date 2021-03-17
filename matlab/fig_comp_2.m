clear;
close all;

load comp_2


hfig = figure(1);
for i=1:length(h0_vec)
    plot(t_vec{i}, J_vec{i})
    hold on;
end
set(gca,'yscale','log');
set(gca,'xscale','log');
set(gca,'FontSize',16);

figure(2);
for i=1:length(h0_vec)
    plot(t_vec{i}(2:end), diff(t_vec{i}))
    hold on;
end
set(gca,'yscale','log');
set(gca,'xscale','log');
set(gca,'FontSize',16);

figure(1)
legend('$h_0=0.001$', '$h_0=0.05$', '$h_0=0.01$', '$h_0=0.1$', '$h_0=0.4$',...
    'Location','SouthWest','interpreter','latex')
ylabel('${f}(R)-{f}(R^*)$','interpreter','latex');
xlabel('$t$','interpreter','latex');

figure(2)
xlabel('$t$','interpreter','latex');
ylabel('$h_k$','interpreter','latex');

conv_x=[4.449 28.94; 4.295 29.95; 4.286 23.14; 4.296 9.956; 4.226 7.916];
conv_y=[0.05102 0.002298; 0.01079 0.0004252; 0.005641 0.0003221; 0.001154 0.0002936; 0.0007076 0.0002483];

plot(conv_x',conv_y','k--');

legend('$h_0=0.001$', '$h_0=0.05$', '$h_0=0.01$', '$h_0=0.1$', '$h_0=0.4$',...
    'Location','NorthEast','interpreter','latex')

text(5.5,5e-2,'$\mathcal{O}(t^{-1.66})$','interpreter','latex','fontsize',18);


diff(log10(conv_y'))./diff(log10(conv_x'))




figure(1);
print('comp_2a','-depsc');
figure(2);
print('comp_2b','-depsc');


!cp comp_2?.eps ../figs

