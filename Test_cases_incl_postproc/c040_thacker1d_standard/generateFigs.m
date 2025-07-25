function generateFigs;

thacker1d_planar;
print(figure(1),'-dpng','-r300',['doc/planar1d.png']);
print(figure(2),'-dpng','-r300',['doc/planar1dwaterlevel.png']);
print(figure(3),'-dpng','-r300',['doc/planar1dvelocity.png']);
close all;