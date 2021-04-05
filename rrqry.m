%RRQRY  Rank-revaeling orthogonal-triangular decomposition.
%   [Q,R,Z] = RRQRY(A), where A is m-by-n, produces an m-by-n upper
%   triangular matrix R and an m-by-m unitary matrix Q so that A(:,Z) =
%   Q*R. Based on Pan&Tang's algorithm number 3.
%
%   [Q,R,Z] = RRQRY(A,0) produces the "economy size" decomposition. If m>n,
%   only the first n columns of Q and the first n rows of R are computed.
%   If m<=n, this is the same as A(:,Z) = Q*R.
%
%   See also RRQRX.