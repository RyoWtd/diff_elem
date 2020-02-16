clear;

order=3;

param = load('./out/out_param.dat');
nu=param(1);
ndf=param(2);
ndt=param(3);
ndr = ndt-ndf
%nd2 = ndt-nd

HH = load('./out/HH.dat');

rr=rank(HH)

p=ndf-rank(HH) %メカニズム自由度
q=nu-rank(HH) %自己釣合自由度

[left, sv, right]=svd(HH);
[n1,n2] = size(right);
ww = right(:,n1);
[m1,m2] = size(left)
force = left(:,rr+1);

svv=diag(sv);
save('./octave/sv-m.dat', 'svv', '-ascii');
save('./octave/right-m.dat', 'ww', '-ascii');
save('./octave/left-m.dat', 'force', '-ascii');

H2 = zeros(nu,ndf,ndf);
t2 = load('./out/HH2.dat');
[nn,mm] = size(t2)

for i = 1:nu
for j = 1:ndf
for k = 1:ndf
H2(i,j,k) = t2((j-1)*ndf+k,i);
end
end
end

if (order>=3)
H3 = zeros(nu,ndf,ndf,ndf);
t3 = load('./out/HH3.dat');
[nn,mm] = size(t3);
for i = 1:nu
for j = 1:ndf
for k = 1:ndf
for m = 1:ndf
H3(i,j,k,m) = t3((j-1)*ndf*ndf+(k-1)*ndf+m,i);
end
end
end
end
end

if (order>=4)
H4 = zeros(nu,ndf,ndf,ndf,ndf);
t4 = load('./out/HH4.dat');
for i = 1:nu
for j = 1:ndf
for k = 1:ndf
for m = 1:ndf
for s = 1:ndf
H4(i,j,k,m,s) = t4((j-1)*ndf*ndf*ndf+(k-1)*ndf*ndf+(m-1)*ndf+s,i);
end
end
end
end
end
end

dd = zeros(ndf,ndf);

g2 = zeros(nu,1);
for i = 1:nu
dd(1:ndf,1:ndf) = H2(i,1:ndf,1:ndf);
g2(i) = -ww'*dd*ww;
end

c2 = zeros(q,1);
for i = 1:q
 ff = left(:,m1-(i-1));
 c2(i) = ff'*g2;
end
c2

ww2 = pinv(HH)*g2;
gg2 = HH*ww2;


save('./octave/g2-m.dat', 'g2', '-ascii');
save('./octave/gg2-m.dat', 'gg2', '-ascii');
save('./octave/c2-m.dat', 'c2', '-ascii');

tmp2 = ww'*ww2;

% idx = load('out_idx.dat');

% www = zeros(ndt,1);
% www2 = zeros(ndt,1);
% for i = 1:ndf
% www(idx(i)) = ww(i);
% www2(idx(i)) = ww2(i);
% end

% save('ww-m.dat', 'www', '-ascii');
% save('ww2-m.dat', 'www2', '-ascii');
save('./octave/ww-m.dat', 'ww', '-ascii');
save('./octave/ww2-m.dat', 'ww2', '-ascii');
% 固定自由度追加までの仮対応

if (order>=3) 
 g3 = zeros(nu,1);
 for i = 1:nu
  for j = 1:ndf
   dd(1:ndf,1:ndf) = H3(i,j,1:ndf,1:ndf);
   ddd = ww'*dd*ww;
   g3(i) = g3(i) - ddd*ww(j);
  end
 end
 
 for i = 1:nu
  dd(1:ndf,1:ndf) = H2(i,1:ndf,1:ndf);
  g3(i) = g3(i) - 3*ww'*dd*ww2;
 end
 
 save('./octave/g3-m.dat', 'g3', '-ascii');
 
 ww3 = pinv(HH)*g3;
 
 gg3=HH*ww3;
 save('./octave/gg3-m.dat', 'gg3', '-ascii');
 
%  www3 = zeros(ndt,1);
%  for i = 1:ndf
%   www3(idx(i)) = ww3(i);
%  end
%  save('./octave/ww3-m.dat', 'www3', '-ascii');
save('./octave/ww3-m.dat', 'ww3', '-ascii');
% 固定自由度追加までの仮対応
 
 c3 = zeros(q,1);
 for i = 1:q
  ff = left(:,m1-(i-1));
  c3(i) = ff'*g3;
 end
 c3
 
 save('./octave/c3-m.dat', 'c3', '-ascii');
 
 
 CCC = [HH,g3];
 rank(CCC)
 
 %[left, sv, right]=svd(CCC);
 
 svv=diag(sv);
 save('./octave/sv2-m.dat', 'svv', '-ascii');
 
 ttt=[left(:,rr+1:m2),g2,g3];
 
 save('./octave/ttt-m.dat', 'ttt', '-ascii');
end

# ここからすべて未修正
% if (order>=4)
%  g4 = zeros(nu,1);
%  for i = 1:nu
%   for j = 1:ndf
%    for k = 1:ndf
%     dd(1:ndf,1:ndf) = H4(i,j,k,1:ndf,1:ndf);
%     ddd = ww'*dd*ww;
%     g4(i) = g4(i) - ddd*ww(j)*ww(k);
%    end
%   end
%  end

%  for i = 1:nu
%   for j = 1:ndf
%    dd(1:ndf,1:ndf) = H3(i,j,1:ndf,1:ndf);
%    ddd = ww2'*dd*ww;
%    g4(i) = g4(i) - 6*ddd*ww(j);
%   end
%  end

%  for i = 1:nu
%   dd(1:ndf,1:ndf) = H2(i,1:ndf,1:ndf);
%   g4(i) = g4(i) - 4*ww3'*dd*ww - 3*ww2'*dd*ww2;
%  end

%  save('g4-m.dat', 'g4', '-ascii');

%  ww4 = pinv(HH)*g4;

%  gg4=HH*ww4;
%  save('gg4-m.dat', 'gg4', '-ascii');

%  www4 = zeros(ndt,1);
%  for i = 1:ndf
%   www4(idx(i)) = ww4(i);
%  end
%  save('ww4-m.dat', 'www4', '-ascii');

%  c4 = zeros(q,1);
%  for i = 1:q
%   ff = left(:,m1-(i-1));
%   c4(i) = ff'*g4;
%  end
%  c4

%  save('c4-m.dat', 'c4', '-ascii');
% 
% end

