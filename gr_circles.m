clear all
% gr_circles.m
% Calculate the number of vertices of a hexagonal tessellation
% that lie on subsequent circles centered 
% at a vertex (OHV), or
% at the center (OHC)
% of one hexagon.
% Plot the lattice and the circles.
% Radii of every vertex around a vertex (OHV) are given by sqrt(A003136)
% Radii of every vertex around the center (OHC) are given by sqrt(A202822)
%
% (c) Szymon Lukaszyk
% email: szymon@patent.pl
% licensed under MIT License
% History
% v1: 17.11.2020
% v2: 08.12.2020 errors due to doubled vertices removed (thanks to Joerg Arndt for hinting them)
% v3: 09.12.2020 improved roundoff errors handling
% Currently this is just a script not a MATLAB function.
% Due to roundoff errors the script may (it will) give inaccurate results for large k.

% FURTHER BUGS ARE ALSO PRESSUMED
% ----------------------------------------------------

% kind of circles
center_at_vertex = 1;
center_at_vertex = 0; % <=> center_at_hex_centre

if center_at_vertex
  OHwhatstr = 'ORIGIN AT HEX VERTEX'; % OHV
else
  OHwhatstr = 'ORIGIN AT HEX CENTER'; % OHC
end
disp( OHwhatstr )

% range of an initial unit square grid p with the origin 0,0
% and size [-rng:1:rng]x[-rng:1:rng]
rng = 40; 

% start with a unit square grid p with the origin 0,0
% ----------------------------------------------------
idx = 1;
for i=-rng:rng % x
  for j=-rng:rng % y
    p(idx,1) = i;
    p(idx,2) = j;    
    idx = idx+1;
  end
end
length(p)
%(2*rng+1)^2
%p
disp 'unit square grid created'

% convert unit square grid to A2 triangle grid
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
length(np)
%floor(length(p)/2) + 1
%np
disp 'A2 triangle grid created'

% create hehagonal grid from A2 grid
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
length(np)
% #formula?
%np
disp 'hehagonal grid created'

% rescale the edge size to a=1 (and center the grid for OHC) 
% ----------------------------------------------------
for k=1:length(np)
  np(k,1) = np(k,1)/sqrt(3);
  np(k,2) = np(k,2)/sqrt(3);
  if ~center_at_vertex
    np(k,1) = np(k,1) + 1;
    np(k,2) = np(k,2);% + sqrt(3)/2;
  end 
end
disp 'hexagonal grid rescaled'

% crop hexagonal grid into a square
% ----------------------------------------------------
nplim = min(max(np))
idx = 1;
for k=1:length(np)
  if ( abs( np(k,1) ) < nplim ) && ( abs( np(k,2) ) < nplim )
    nnp(idx,1) = np(k,1);
    nnp(idx,2) = np(k,2);   
    idx = idx+1;  
  end
end
np = nnp;
clear nnp
disp 'hexagonal grid cropped'

hx=np;

% calculate radii rn between (0,0) and every hx
% ----------------------------------------------------
for k=1:length(hx)
  rn(k) = ( hx(k,1)^2 + hx(k,2)^2 )^(1/2);
end
rn = sort(rn);
[aa, bb, rn]= find( rn ); % get rid of zeros for the case of OHV
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
disp 'out of the grid radii removed'

% calculate differences between radii
% ----------------------------------------------------
for k=2:length(rn)
  diff(k-1) = rn(k)-rn(k-1);
end

diff=sort(diff);
for k=1:length(diff)
  if diff(k) > 10^-13 % set with rng=300 %!! BUT THIS VALUE IS SURELY DERIVABLE FROM diff !!
     kthr1 = k-1;  % max diff due to numerical errors
     kthr2 = k;    % min diff due to graphene structure
     break
  end
end
disp 'differences calculated'

diff(kthr1), diff(kthr2)
%[aa, bb, diff]= find( diff ); % get rid of zeros
%figure, plot( diff(1:length(diff)-5) ), title(OHwhatstr);

% collect radii rnd having the same length
% ----------------------------------------------------
rnd(1) = rn(1); % rn(1) = 1
idx = 2;

%rnd=0;
%idx = 1;
for k=1:length(rn)
  % take care of roundoff errors
  findres = find( abs( rnd-rn(k) ) < diff(kthr2) );
  
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
for k=1:length(rn)
  if abs( rnd(idx)-rn(k) ) < diff(kthr2)
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
disp 'vertices counted'

% gather unique count values
% ----------------------------------------------------
countun(1) = count(1);
idx = 2;
for k=2:length(count)
  findres = find( countun == count(k) );
  if isempty( findres )
    countun(idx) = count(k);
    idx = idx+1;         
  end
end
countun=sort(countun);
countun
disp 'unique count values gathered'

% validity check 1
% ----------------------------------------------------
sum(count) == length(rn)

% validity check 2
for k=1:length(count)
  if mod( count(k), 3 ) ~= 0
    disp( sprintf('Error for %d circle yielding %d count', k, count(k) ) )
  end
end   

%return

% draw circles around a vertex 0,0
% ----------------------------------------------------
figure
scatter(hx(:,1), hx(:,2),'.')
hold on

colkind = 'r'; 
%colkind = 'b'; 
for k=1:length(rnd)
  th = 0:pi/50:2*pi;
  xunit = rnd(k) * cos(th);
  yunit = rnd(k) * sin(th);
  plot(xunit, yunit, colkind);
  
  if colkind == 'r'
    colkind = 'b';
    text( rnd(k), 0.2, sprintf('%d', k) )    
  else 
    colkind = 'r';  
    text( 0.2, rnd(k), sprintf('%d', k) )        
    
  end
end

return
