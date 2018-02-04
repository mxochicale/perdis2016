function [pqr, pqd, pqi] = kd_search(y,r,tree,yp)
% KD_SEARCH     search kd_tree for r nearest neighbors of point yq.
%               Need to partition original data using kd_partition.m
%               [r,tree] = kd_partition(ym); ym--original data
%
% USE: [pqr, pqd, pqi] = kd_search(y,r,tree,yp);
%
% INPUT:
%       y: array of query points in columnwise form (size(y,1)=d)
%       r: requested number of nearest neighbors to the query point yq
%       tree: kd_tree constructed using kd_partition.m
%       yp: partitioned (ordered) set of data that needs to be searched
%           (using my kd_partirion you want to input ym(:,indx), where ym 
%           is the data used for partitioning and indx sorted index of ym)
%
% OUTPUT:
%       pqr: cell array of the r nearest neighbors of y in yp 
%       pqd: cell array of the corresponding distances
%       pqi: cell array of the indices of r nearest neighbors of y in yp to
%            get indexes for y, use indx(pqi{1})
%
% copyright (c) and written by David Chelidze, February 02, 2009. updated
% 11/26/2012
 
% check inputs
if nargin < 4
    error('Need four input variables to work')
end
 
% declare global variables for all subfunctions
global yq qri qr qrd b_lower b_upper
 
n = size(y,2);
pqr = cell(n,1);
pqd = cell(n,1);
pqi = cell(n,1);
 
for k = 1:n,
    yq = y(:,k);
    qrd = []; % initialize array for r distances
    qr = []; % initialize r nearest neighbor points
    qri = []; % initialize index of r nearest neighbors
 
    % set up the box bounds, which start at +/- infinity (whole k-d space)
    b_upper = Inf*ones(size(yq));
    b_lower = -b_upper;
 
    kdsrch(1,r,tree,yp); % start the search from the first node
    pqr{k} = qr;
    pqd{k} = sqrt(qrd);
    pqi{k} = qri;
end
 
%--------------------------------------------------------------------------            
function kdsrch(node,r,tree,yp)
% KDSRCH    actual k-d search
% this drills down the tree to the end node and updates the 
% nearest neighbors list with new points
%
% INPUT: starting node number, and kd_partition data
%
% copyright (c) and written by David Chelidze, February 02, 2009, 
% modified 11/26/2012
 
global yq qri qr qrd b_lower b_upper
 
if tree.left(node) == 0 % this is a terminal node: update the nns
 
    qri = [qri tree.first(node):tree.last(node)]; % update nns indexes
    qr = yp(:,qri); % current list of nns including points in this bin
    qrd = sum((qr - yq*ones(1,size(qr,2))).^2,1); % distances squared
    [qrd, indx] = sort(qrd); % sorted distances squared and their index
    qr = qr(:,indx); % sorted list of nearest neighbors
    qri = qri(indx); % sorted list of indexes
    if size(qr,2) > r % truncate to the first r points
        qrd = qrd(1:r);
        qr = qr(:,1:r);
        qri = qri(1:r);
    end
    % be done if all points are with this box
    if within(yq, b_lower, b_upper, qrd(end)) && numel(qri) == r
        return
    end % otherwise (during backtracking) WITHIN will always return 0.   
 
else
 
    d = tree.splitdim(node); % split dimension for current node
    p = tree.splitval(node); % the corresponding split value
 
    % first determine which child node to search
    if yq(d) <= p % need to search left child node
        tmp = b_upper(d);
        b_upper(d) = p;
        kdsrch(tree.left(node),r,tree,yp);
        b_upper(d) = tmp;
    else % need to search the right child node
        tmp = b_lower(d);
        b_lower(d) = p;
        kdsrch(tree.right(node),r,tree,yp);
        b_lower(d) = tmp;
    end
 
    % check if other nodes need to be searched (backtracking)
    if yq(d) <= p
        tmp = b_lower(d);
        b_lower(d) = p;
        if overlap(yq, b_lower, b_upper, qrd(end)) 
            % need to search the right child node
            kdsrch(tree.right(node),r,tree,yp);
        end
        b_lower(d) = tmp;
    else
        tmp = b_upper(d);
        b_upper(d) = p;
        if overlap(yq, b_lower, b_upper, qrd(end)) 
            % need to search the left child node
            kdsrch(tree.left(node),r,tree,yp);
        end
        b_upper(d) = tmp;
    end % when all the other nodes are searched
    
end
 
% see if we should terminate search
if within(yq, b_lower, b_upper, qrd(end)) && numel(qri) == r
    return
end % otherwise (during backtracking) WITHIN will always return 0.
 
 
%--------------------------------------------------------------------------            
function flag = within(yq, b_lower, b_upper, ball)
% WITHIN    check if additional nodes need to be searched (i.e. if the ball
%  centered at yq and containing all current nearest neighbors overlaps the
%  boundary of the leaf box containing yq)
%
% INPUT:
%   yq: query point
%   b_lower, b_upper: lower and upper bounds on the leaf box
%   ball: square of the radius of the ball centered at yq and containing
%         all current r nearest neighbors
% OUTPUT:
%   flag: 1 if ball does not intersect the box, 0 if it does
%
% Modified by David Chelidze on 02/03/2009.
 
if ball <= min([abs(yq-b_lower)', abs(yq-b_upper)'])^2
    % ball containing all the current nn is inside the leaf box (finish)
    flag = 1;
else % ball overlaps other leaf boxes (continue recursive search)
    flag = 0; 
end
 
%--------------------------------------------------------------------------            
function flag = overlap(yq, b_lower, b_upper, ball)
% OVERLAP   check if the current box overlaps with the ball centered at yq
%   and containing all current r nearest neighbors.
%
% INPUT:
%   yq: query point
%   b_lower, b_upper: lower and upper bounds on the current box
%   ball: square of the radius of the ball centered at yq and containing
%         all current r nearest neighbors
% OUTPUT:
%   flag: 0 if ball does not overlap the box, 1 if it does
%
% Modified by David Chelidze on 02/03/2009.
 
il = find(yq < b_lower); % index of yq coordinates that are lower the box 
iu = find(yq > b_upper); % index of yq coordinates that are upper the box
% distance squared from yq to the edge of the box
dst = sum((yq(il)-b_lower(il)).^2,1)+sum((yq(iu)-b_upper(iu)).^2,1);
if dst >= ball % there is no overlap (finish this search)    
    flag = 0;
else % there is overlap and the box needs to be searched for nn
    flag = 1;
end


