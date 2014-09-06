function Q = orth2(A,tol)
if isempty(A)
    Q=A;
    return;
end
[U,S] = svd(A,0);
[m,n] = size(A);
if m > 1, s = diag(S);
   elseif m == 1, s = S(1);
   else s = 0;
end
r = sum(s > max(m,n) * max(s(1),1) * tol);
Q = U(:,1:r);