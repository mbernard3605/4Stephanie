function [lam, lds, pcs, per] = eof_routine(data, nkp);
%
%  [lam, lds, pcs, per] = eof_routine(data);
%
%  This function assumes that [ntim, npts] = size(data);
%  that is, each column is a time series for a different point.
%
if nargin == 1;
  nkp = 10;
end

if ndims(data) ~= 2;
  szdat = size(data);
  data = reshape(data, szdat(1), prod(szdat(2:end)));
end

[m, n] = size(data);

rot = (m < n);
if rot;
  data = data'; %[m, n] = size(data);
end
%rot = 0;
%[m,n]=size(data);
c = data' * data / (m-1);
[lds, lam] = eig(c);
l = diag(lam);

[lam, k] = sort(l'); lds = lds(:,k);

lam = fliplr(lam);
lds = fliplr(lds); lds = lds(:,1:nkp);
per = lam * 100 / sum(lam);
pcs = data*lds;

if rot;
  tem = pcs;
  pcs = lds;
  lds = tem;
  for i = 1:nkp;
    pcs(:,i) = pcs(:,i) / sqrt(lam(i)*(m-1));
    lds(:,i) = lds(:,i) * sqrt(lam(i)*(m-1));
  end
end
