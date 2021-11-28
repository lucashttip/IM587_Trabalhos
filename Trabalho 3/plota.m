function plota(x,y,labelx, labely, filename)


  f = figure;
  set (f,'defaultaxesfontsize', 13);
%  set (f,'position',[300 150 800 600]);
  
  ax = axes;
  
  plot(x,y,"linewidth",2);
  set(ax, "linewidth", 1.5, "labelfontsizemultiplier",1.2);
  ylabel(ax,labely);
  xlabel(ax, labelx);

  
%  set (f,'position',[300 150 800 600]);
  print (f, filename);
##  print -dpng -color plot1.png

end