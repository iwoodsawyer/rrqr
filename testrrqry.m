n = 42;
m = 41;

A = randn(m,n);

tic
[Q,R,p,r] = rrqry(A);
toc
norm(A(:,p)-Q*R)

tic
[Q,R,p,r] = rrqry(A,0);
toc
norm(A(:,p)-Q*R)

n = 42;
m = 61;

A = randn(m,n);

tic
[Q,R,p,r] = rrqry(A);
toc
norm(A(:,p)-Q*R)

tic
[Q,R,p,r] = rrqry(A,0);
toc
norm(A(:,p)-Q*R)

n = 61;
m = 42;

A = randn(m,n);

tic
[Q,R,p,r] = rrqry(A);
toc
norm(A(:,p)-Q*R)

tic
[Q,R,p,r] = rrqry(A,0);
toc
norm(A(:,p)-Q*R)