function generateFigs;

% Beware! The run time shortened for the testbench
%         A fully developed outcome requires about 3000 hours of simulated
%         time.

% Print
poiseuillepartialslip
print(figure(1),'-dpng','-r300','doc/poiseuillepartialslip.png');
close all;