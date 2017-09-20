function c = color_interp(q,cols,range)
  % colormap is cols,
  % q is between 0 and 1

  x = linspace(range(1),range(2),size(cols,1));
c1 = interp1(x,cols(:,1),q);
c2 = interp1(x,cols(:,2),q);
c3 = interp1(x,cols(:,3),q);

c = [c1 c2 c3];
