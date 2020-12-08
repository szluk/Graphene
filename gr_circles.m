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
% v1: 17.11.2020
% v2: 08.12.2020 errors due to double vertices removed (thanks to Joerg Arndt)
% Currently this is just a script not a MATLAB function.
% Due to roundoff errors it may (it will) give inaccurate results for large k.

% FURTHER SERIOUS BUGS ARE PRESSUMED

% kind of circles
center_at_vertex = 1;
center_at_vertex = 0; % <=> center_at_ring_centre

if center_at_vertex
  disp 'ORIGIN AT HEX VERTEX'
else
  disp 'ORIGIN AT HEX CENTER'
end

% range of an initial unit square grid p with the origin 0,0
% and size [-rng:1:rng]x[-rng:1:rng]
rng = 100; 
rng = 40; 

% start with a unit square grid p with the origin 0,0
% ----------------------------------------------------
idx = 1;
for i=-rng:rng % x
  for j=-rng:rng % y
    %sprintf('%d,%d', i, j)
    p(idx,1) = i;
    p(idx,2) = j;    
    idx = idx+1;
  end
end
length(p) % 
%(2*rng+1)^2
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
length(np) % |np| = floor(length(p)/2) + 1
%floor(length(p)/2) + 1
%np
disp 'A2 triangle grid np created'

% create hehagon grid from A2 grid
% ----------------------------------------------------
rngp_flg=1;
idx = 1;
k=1;
while k<=length(np)
  for l=0:2*rng
    if k+l>length(np)
      break
    end
    nnp(idx,:) = np(k+l,:);
    %sprintf('%d,%d', k+l, idx)
    idx = idx+1;
  end
  if rngp_flg
    k=k+l+rng+2;
    rngp_flg=0;
  else
    k=k+l+rng+1;
    rngp_flg=1;  
  end  
end
np = nnp;
clear nnp
disp 'hexagon grid np created'

% rescale to a=1 and center
% ----------------------------------------------------
for k=1:length(np)
  np(k,1) = np(k,1)/sqrt(3);
  np(k,2) = np(k,2)/sqrt(3);
  if ~center_at_vertex
    np(k,1) = np(k,1) + 1;
    np(k,2) = np(k,2);% + sqrt(3)/2;
  end 
end
disp 'hexagon grid np rescaled'

% crop A2 triangle grid np into a square
% ----------------------------------------------------
nplim = min(max(np))
idx = 1;
for k=1:length(np)
  if ( abs( np(k,1) ) < nplim ) && ( abs( np(k,2) ) < nplim )
    nnp(idx,1) = np(k,1);
    nnp(idx,2) = np(k,2);   
    %sprintf('%d,(%d,%d)', nnp(idx,1), nnp(idx,2))
    idx = idx+1;  
  end
end
np = nnp;
clear nnp
disp 'A2 triangle grid cropped'

figure
scatter(np(:,1), np(:,2),'.')
hold on
%return

hx=np;

% calculate radii rn between (0,0) and every hx
% ----------------------------------------------------
for k=1:length(hx)
  rn(k) = ( hx(k,1)^2 + hx(k,2)^2 )^(1/2);
end
rn = sort(rn);
%rn
disp 'radii calculated'

% remove radii rn that are out of the grid
% ----------------------------------------------------
idx = 1;
for k=1:length(rn)
  if rn(k) <= min(abs(min(hx))-.9) % out of the grid with a safe margin  
    rn1(idx) = rn(k);
    idx = idx+1;         
  end
end
rn = rn1;
clear rn1
disp 'radii removed'

% collect radii rnd having the same length
% ----------------------------------------------------
rnd = 0;
idx = 1;
rncnt = 0;
for k=1:length(rn)
  rncnt = rncnt+1;
  
  % take care of roundoff errors
  %findres = find( abs( rnd-rn(k) ) < 10^(-13)  ); % 10^(-13) is OK up to rng=250, k=9999 (A003136 table)
  %findres = find( abs( rnd-rn(k) ) < 10^(-12)  ); % 10^(-12) is OK to A202822 table
  findres = find( abs( rnd-rn(k) ) < 10^(-13)  );
  if isempty( findres )
    rnd(idx) = rn(k);
    idx = idx+1;         
  end
end
rnd = sort(rnd);
disp 'same radii created'

% count vertices at the same radius rn
% ----------------------------------------------------
idx = 1;
count(idx) = 0;

if center_at_vertex
  kmin=2
else
  kmin=1
end

for k=kmin:length(rn)
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

count

% validity check
sum(count) == length(rn)

% draw circles around a vertex 0,0
% ----------------------------------------------------
colkind = 'r'; 
%colkind = 'b'; 
for k=1:length(rnd)
  th = 0:pi/50:2*pi;
  xunit = rnd(k) * cos(th);
  yunit = rnd(k) * sin(th);
  plot(xunit, yunit, colkind);
  %text( rnd(k), 0, sprintf('\\leftarrow %d', k) )      

  
  if colkind == 'r'
    colkind = 'b';
    text( rnd(k), 0.2, sprintf('%d', k) )    
  else 
    colkind = 'r';  
    text( 0.2, rnd(k), sprintf('%d', k) )        
    
  end
end

return
