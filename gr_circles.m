clear all
% gr_circles.m
% Calculate the number of vertices in subsequent circles centered 
% at a vertex or at the centre of a ring of a hexagonal tessellation.
% Radii rc of every vertex around a vertex are given by sqrt(A003136)
% Radii rc of every vertex around a ring center are given by sqrt(A202822)
% ----------------------------------------------------

% Plot the lattice and the circles.
%
% (c) Szymon Lukaszyk
% email: szymon@patent.pl
% licensed under MIT License
% History
% v1:0 17.11.2020
% Currently this is just a script not a MATLAB function.
% Due to roundoff errors it may (it will) give inaccurate results for large k.

% kind of circles
% center_at_vertex = 1;
center_at_vertex = 0; % <=> center_at_ring_centre

% range of an initial unit square grid p with the origin 0,0
% and size [-rng:1:rng]x[-rng:1:rng]
rng = 100; 

% start with a unit square grid p with the origin 0,0
% ----------------------------------------------------
idx = 1;
for i=-rng:rng % x
  for j=-rng:rng % y
    p(idx,:) = [i j];
    idx = idx+1;
  end
end
%p
disp 'unit square grid p created'

% convert unit square grid to A2 triangle grid np
% ----------------------------------------------------
x_dlt = sqrt(3)/2;
y_dlt = 3/2;

idx = 1;
for k=1:length(p)
    if rem( abs(p(k,1)) + abs(p(k,2)), 2 ) == 0 % even check
      np(idx,1) = p(k,1)*x_dlt;  
      np(idx,2) = p(k,2)*y_dlt;      
      idx = idx+1;
    end
end
%np
disp 'A2 triangle grid np created'

% crop A2 triangle grid np into a square
% ----------------------------------------------------
nplim = min(max(np));
idx = 1;
for k=1:length(np)
  if abs( np(k,2) ) < nplim
    nnp(idx,:) = np(k,:);
    idx = idx+1;  
  end
end
np = nnp;
disp 'A2 triangle grid cropped'

% create hehagon around every vertes of the A2 grid
% and remove duplicates
% ----------------------------------------------------
%              0, 1 (2)
%          /         \
%-p3/2, 1/2 (1)       p3/2, 1/2 (3)
%         |    x,y    |
%-p3/2,-1/2 (6)       p3/2,-1/2 (4)
%          \         /
%              0,-1 (5)
hx = [0 0];
idx = 1;
for k=1:length(np)
  x = np(k,1);
  y = np(k,2);  
  for i=1:6
    switch i
      case 1
        xx = x-sqrt(3)/2; 
        yy = y+1/2;
      case 2
        xx = x; 
        yy = y+1;
      case 3
        xx = x+sqrt(3)/2; 
        yy = y+1/2;
      case 4
        xx = x+sqrt(3)/2; 
        yy = y-1/2;
      case 5
        xx = x; 
        yy = y-1;
      case 6
        xx = x-sqrt(3)/2; 
        yy = y-1/2;
      end
      % take care of roundoff errors
      if isempty( find( abs(hx(:,1)-xx)<10^(-15) & abs(hx(:,2)-yy)<10^(-14) ) ) % single & here      
         idx = idx+1;         
         hx(idx,1) = xx;
         hx(idx,2) = yy;         
      end
  end
end
disp 'hex created'

figure
hold on
grid on
%scatter(np(:,1), np(:,2),'o')

if center_at_vertex 
  % create new heh grid with the origin in the vertex
  % ----------------------------------------------------
  for k=2:length(hx) % hx(0) is the origin (0,0)
    x = hx(k,1) + sqrt(3)/2;
    y = hx(k,2) + 1/2;  
    hxn(k-1, 1) = x;
    hxn(k-1, 2) = y;  
  end
  hx = hxn;
  disp 'ORIGIN AT HEX VERTEX'  
else
  disp 'ORIGIN AT HEX CENTER'  
end  
%hx
scatter(hx(:,1), hx(:,2),'.r')

% populate all radii first
% ----------------------------------------------------
for k=1:length(hx)
  rn(k) = ( hx(k,1)^2 + hx(k,2)^2 )^(1/2);
end
rn = sort(rn);
%rn
disp 'radii populated'

% truncate radii out of the grid
% ----------------------------------------------------
idx = 1;
for k=1:length(rn)
  if rn(k) <= min(abs(min(hx))-.9) % out of the grid with a safe margin  
    rn1(idx) = rn(k);
    idx = idx+1;         
  end
end
rn = rn1;
disp 'radii truncated'

% collect the same radii
% ----------------------------------------------------
rnd = 0;
idx = 1;
rncnt = 0;
for k=1:length(rn)
  rncnt = rncnt+1;
  
  % take care of roundoff errors
  %findres = find( abs( rnd-rn(k) ) < 10^(-13)  ); % 10^(-13) is OK up to rng=250, k=9999 (A003136 table)
  findres = find( abs( rnd-rn(k) ) < 10^(-12)  ); % 10^(-12) is OK to A202822 table
  if isempty( findres )
    rnd(idx) = rn(k);
    idx = idx+1;         
  end
end
rnd = sort(rnd);
disp 'same radii created'

% count vertices at the same radius
% ----------------------------------------------------
idx = 1;
count(idx) = 0;
for k=2:length(rn) % rn(1) = 0
  if abs( rnd(idx)-rn(k) ) < 10^(-9)   
    count(idx) = count(idx)+1;
  else
    idx = idx+1;
    if idx > length(rnd)
      break
    end
    count(idx) = 1;
  end
end

cntsm(1) = count(1);
for k=2:length(count)
  cntsm(k) = cntsm(k-1) + count(k);
end

count

% validity check
sum(count) == length(rn)-1

% draw circles around a vertex 0,0
% ----------------------------------------------------
colkind = 'r'; 
%colkind = 'b'; 
for k=1:length(rnd)
  th = 0:pi/50:2*pi;
  xunit = rnd(k) * cos(th);
  yunit = rnd(k) * sin(th);
  plot(xunit, yunit, colkind);
  if colkind == 'r'
    colkind = 'b';
  else 
    colkind = 'r';  
  end
end

return
