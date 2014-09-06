function Z = null2(A,tol)

[m,n] = size(A);
[U,S,V]=svd(A,0);
if m > 1
    s = diag(S);
  elseif m == 1
      s = S(1);
  else s = 0;
end
r = sum(s > max(m,n) * max(s(1),1) * tol);
Z = V(:,r+1:n);
end
