function [fnn,Dx,Dz] = false_neighbors(x, m, D, n, R, s) 
% FALSE_NEIGHBORS estimate appropriate embedding dimension by estimating 
% false nearest neighbors in a reconstructed phase space based on a 
% delay coordinate embedding similar to: 
% "Determining embedding dimension for phase-space 
% reconstruction using a geometrical construction", 
% M. B. Kennel, R. Brown, and H.D.I. Abarbanel, 
% Physical Review A, Vol 45, No 6, 15 March 1992, 
% pp 3403-3411. 
% the method was further improved to account for (d+1) nn components 
% that were larger than mean random distance, nn distances that were 
% larger than mean random distances. 
% for more details, please see this post: 
% http://egr.uri.edu/nld/2013/12/16/notes-on-false-nearest-neighbors/ 
% 
% CALL: 
% [fnn,Dx,Dz] = false_neighbors(x, m, D, n, r, t); 
% 
% INPUT: % x - scalar time series 
% m - delay time 
% D - maximum embedding dimension to consider 
% n - number of points to consider 
% R - (optional) is used to detect fnn when z(k+1) is greater than R*d(k) 
% the default value is 5 and KBA suggests (15-20) 
% s - (between 0 and 1, optional) determines the threshold for z(k+1) 
% beeing greater than std of uncorr random noise (default is 0.5) 
% 
% OUTPUT: 
% fnn - ratio of false nearest neighbors in each dimension 
% Dx - cell array of distances to nearest neighbors in each dimension 
% Dz - cell array of distances between d+1 coordinates of nns in d-dims 
% 
% copyright 2013 David Chelidze 

if nargin < 6 || isempty(s) % use default threshold 
    s = 0.5; 
end 
if nargin < 5 || isempty(R) % use default threshold 
    R = 5; 
end 
if nargin < 4 % need time series error('ERROR: need at least four input variables') 
elseif min(size(x)) == 1 % normalize input time series 
    x = (x(:)' - mean(x))/std(x); 
    nm = size(x, 2) - D*m; 
    n = min(nm, n); 
else % cannot deal with vector-valued time series 
    error('ERROR: x needs to be a scalar time series') 
end 

% reconstruct the phase space 
y = zeros(D+1, n); 
for i = 0:D y(i+1,:) = x((1:n) + i*m); 
end
fnn = zeros(1,D); % initialize fnn array 
Dx= cell(1,D); 
Dz = Dx; 
for d = 1:D 
    % find the false nearest neighbors for each embedding dimension 
    fprintf('Estimating false nearest neighbors for d = %d\n', d) 
    % partition all phase space points for fast nn searching 
    [tree, r] = kd_partition(y(1:d, :), 256); 
    % find the nearest neighbor for each of the phase space point 
    [~, nnd, nni] = kd_search(y(1:d, :), 2, tree, y(1:d, r)); 
    nni = cell2mat(nni); % convert to double 
    nni = r(nni(:,2)); % indexes of nn point 
    nnd = cell2mat(nnd); % convert to double 
    nnd = nnd(:,2); % distances to nn point 
    dz = abs(y(d+1, :) - y(d+1, nni))'; % delta lengths in (d+1)-dimension 
    if nargout > 1 % stode dnns 
        Dx{d} = nnd; % for randn data has d-d.o.f. chi distribution 
elseif nargout == 2 % store dzs 
    Dz{d} = dz; % for randn data dz has chi distribution at each dim 
    end
    dr = dz./nnd; % ratio of (d+1)-dim delta lengths over d-dim distances 
    % check if dz is greater than s*std of data or if the dr > R 
    fnn(d) = sum(dz > s | dr > R)/length(dz); % fnns in d-dimensions 
end
if nargout == 0 % plot the results 
    plot(1:D,100*fnn,'o-') 
    axis([0 D 0 100]) 
    set(gca,'fontsize',12) 
    xlabel('Embedding Dimension (size)') 
    ylabel('FNNs (%)')
end



%%-------------------------------------------------------------------------
 function [pqr, pqd, pqi] = kd_search(y,r,tree,yp) 
     % KD_SEARCH search kd_tree 
for r nearest neighbors of point yq. % Need to partition original data using 
kd_partition.m % [r,tree] = kd_partition(ym); ym--original data % % USE: [pqr, 
pqd, pqi] = kd_search(y,r,tree,yp); % % INPUT: % y: array of query points in 
columnwise form (size(y,1)=d) % r: requested number of nearest neighbors to the 
query point yq % tree: kd_tree constructed using kd_partition.m % yp: 
partitioned (ordered) set of data that needs to be searched % (using my 
kd_partirion you want to input ym(:,indx), where ym % is the data used for 
partitioning and indx sorted index of ym) % % OUTPUT: % pqr: cell array of the r 
nearest neighbors of y in yp % pqd: cell array of the corresponding distances % 
pqi: cell array of the indices of r nearest neighbors of y in yp to % get 
indexes for y, use indx(pqi{1}) % % copyright (c) and written by David Chelidze, 
February 02, 2009. updated % 11/26/2012 % check inputs if nargin < 4 error('Need 
four input variables to work') end % declare global variables for all 
subfunctions global yq qri qr qrd b_lower b_upper n = size(y,2); pqr = 
cell(n,1); pqd = cell(n,1); pqi = cell(n,1); for k = 1:n, yq = y(:,k); qrd = []; 
% initialize array for r distances qr = []; % initialize r nearest neighbor 
points qri = []; % initialize index of r nearest neighbors % set up the box 
bounds, which start at +/- infinity (whole k-d space) b_upper = 
Inf*ones(size(yq)); b_lower = -b_upper; kdsrch(1,r,tree,yp); % start the search 
from the first node pqr{k} = qr; pqd{k} = sqrt(qrd); pqi{k} = qri; end
 %--------------------------------------------------------------------------
 function kdsrch(node,r,tree,yp) % KDSRCH actual k-d search % this drills down 
the tree to the end node and updates the % nearest neighbors list with new 
points % % INPUT: starting node number, and kd_partition data % % copyright (c) 
and written by David Chelidze, February 02, 2009, % modified 11/26/2012 global 
yq qri qr qrd b_lower b_upper if tree.left(node) == 0 % this is a terminal node: 
update the nns qri = [qri tree.first(node):tree.last(node)]; % update nns 
indexes qr = yp(:,qri); % current list of nns including points in this bin qrd = 
sum((qr - yq*ones(1,size(qr,2))).^2,1); % distances squared [qrd, indx] = 
sort(qrd); % sorted distances squared and their index qr = qr(:,indx); % sorted 
list of nearest neighbors qri = qri(indx); % sorted list of indexes if 
size(qr,2) > r % truncate to the first r points qrd = qrd(1:r); qr = qr(:,1:r); 
qri = qri(1:r); end % be done if all points are with this box if within(yq, 
b_lower, b_upper, qrd(end)) && numel(qri) == r return end % otherwise (during 
backtracking) WITHIN will always return 0. else d = tree.splitdim(node); % split 
dimension for current node p = tree.splitval(node); % the corresponding split 
value % first determine which child node to search if yq(d) <= p % need to 
search left child node tmp = b_upper(d); b_upper(d) = p; 
kdsrch(tree.left(node),r,tree,yp); b_upper(d) = tmp; else % need to search the 
right child node tmp = b_lower(d); b_lower(d) = p; 
kdsrch(tree.right(node),r,tree,yp); b_lower(d) = tmp; end % check if other nodes 
need to be searched (backtracking) if yq(d) <= p tmp = b_lower(d); b_lower(d) = 
p; if overlap(yq, b_lower, b_upper, qrd(end)) % need to search the right child 
node kdsrch(tree.right(node),r,tree,yp); end b_lower(d) = tmp; else tmp = 
b_upper(d); b_upper(d) = p; if overlap(yq, b_lower, b_upper, qrd(end)) % need to 
search the left child node kdsrch(tree.left(node),r,tree,yp); end b_upper(d) = 
tmp; end % when all the other nodes are searched end % see if we should 
terminate search if within(yq, b_lower, b_upper, qrd(end)) && numel(qri) == r 
return end % otherwise (during backtracking) WITHIN will always return 0.
 %-------------------------------------------------------------------------- 
function flag = within(yq, b_lower, b_upper, ball) % WITHIN check if additional 
nodes need to be searched (i.e. if the ball % centered at yq and containing all 
current nearest neighbors overlaps the % boundary of the leaf box containing yq)
 % % INPUT: % yq: query point % b_lower, b_upper: lower and upper bounds on the 
leaf box % ball: square of the radius of the ball centered at yq and containing
 % all current r nearest neighbors % OUTPUT: % flag: 1 if ball does not 
intersect the box, 0 if it does % % Modified by David Chelidze on 02/03/2009. if 
ball <= min([abs(yq-b_lower)', abs(yq-b_upper)'])^2 % ball containing all the 
current nn is inside the leaf box (finish) flag = 1; else % ball overlaps other 
leaf boxes (continue recursive search) flag = 0; end
 %--------------------------------------------------------------------------
 function flag = overlap(yq, b_lower, b_upper, ball) % OVERLAP check if the 
current box overlaps with the ball centered at yq % and containing all current r 
nearest neighbors. % % INPUT: % yq: query point % b_lower, b_upper: lower and 
upper bounds on the current box % ball: square of the radius of the ball 
centered at yq and containing % all current r nearest neighbors % OUTPUT: % 
flag: 0 if ball does not overlap the box, 1 if it does % % Modified by David 
Chelidze on 02/03/2009. il = find(yq < b_lower); % index of yq coordinates that 
are lower the box iu = find(yq > b_upper); % index of yq coordinates that are 
upper the box % distance squared from yq to the edge of the box dst = 
sum((yq(il)-b_lower(il)).^2,1)+sum((yq(iu)-b_upper(iu)).^2,1); if dst >= ball % 
there is no overlap (finish this search) flag = 0; else % there is overlap and 
the box needs to be searched for nn flag = 1; end
 %--------------------------------------------------------------------------
 function [kd_tree, r] = kd_partition(y, b, c) % KD_PARTITION Create a kd-tree 
and partitioned database for efficiently % finding the nearest neighbors to a 
point in a % d-dimensional space. % % USE: [kd_tree, r] = kd_partition(y, b, c);
 % % INPUT: % y: original multivariate data (points arranged columnwise, 
size(y,1)=d) % b: maximum number of distinct points for each bin (default is 
100) % c: minimum range of point cluster within a final leaf (default is 0) % % 
OUTPUT: % kd_tree structure with the following variables: % splitdim: dimension 
used in splitting the node % splitval: corresponding cutting point % first & 
last: indices of points in the node % left & right: node #s of consequent 
branches to the current node % r: sorted index of points in the original y 
corresponding to the leafs % % to find k-nearest neighbors use kd_search.m % % 
copyrighted (c) and written by David Chelidze, January 28, 2009. % check the 
inputs if nargin == 0 error('Need to input at least the data to partition')
 elseif nargin > 3 error('Too many inputs') end % initializes default variables 
if needed if nargin < 2 b = 100; end if nargin < 3 c = 0; end [d, n] = size(y); 
% get the dimension and the number of points in y r = 1:n; % initialize original 
index of points in y % initializes variables for the number of nodes and the 
last node node = 1; last = 1; % initializes the first node's cut dimension and 
value in the kd_tree kd_tree.splitdim = 0; kd_tree.splitval = 0; % initializes 
the bounds on the index of all points kd_tree.first = 1; kd_tree.last = n; % 
initializes location of consequent branches in the kd_tree kd_tree.left = 0;
 kd_tree.right = 0; while node <= last % do until the tree is complete % specify 
the index of all the points that are partitioned in this node segment = 
kd_tree.first(node):kd_tree.last(node); % determines range of data in each 
dimension and sorts it [rng, index] = sort(range(y(:,segment),2)); % now 
determine if this segment needs splitting (cutting) if rng(d) > c && 
length(segment)>= b % then split yt = y(:,segment); rt = r(segment); [sorted, 
sorted_index] = sort(yt(index(d),:)); % estimate where the cut should go lng = 
size(yt,2); cut = (sorted(ceil(lng/2)) + sorted(floor(lng/2+1)))/2; L = (sorted 
<= cut); % points to the left of cut if sum(L) == lng % right node is empty L = 
(sorted < cut); % decrease points on the left cut = (cut + max(sorted(L)))/2; % 
adjust the cut end % adjust the order of the data in this node y(:,segment) = 
yt(:,sorted_index); r(segment) = rt(sorted_index); % assign appropriate split 
dimension and split value kd_tree.splitdim(node) = index(d); 
kd_tree.splitval(node) = cut; % assign the location of consequent bins and 
kd_tree.left(node) = last + 1; kd_tree.right(node) = last + 2; % specify which 
is the last node at this moment last = last + 2; % initialize next nodes cut 
dimension and value in the kd_tree % assuming they are terminal at this point 
kd_tree.splitdim = [kd_tree.splitdim 0 0]; kd_tree.splitval = [kd_tree.splitval 
0 0]; % initialize the bounds on the index of the next nodes kd_tree.first = 
[kd_tree.first segment(1) segment(1)+sum(L)]; kd_tree.last = [kd_tree.last 
segment(1)+sum(L)-1 segment(lng)]; % initialize location of consequent branches 
in the kd_tree % assuming they are terminal at this point kd_tree.left = 
[kd_tree.left 0 0]; kd_tree.right = [kd_tree.right 0 0]; end % the splitting 
process % increment the node node = node + 1; end % the partitioning