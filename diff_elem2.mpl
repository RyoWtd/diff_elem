restart;
Digits:=18;

with(LinearAlgebra):

ord:=3: #次数
tyl:=3: #テイラー展開次数

nue:=20: #要素非適合量数
# (12成分+ 8成分 e1_1,e2_1,e3_1,e2_1,e2_2,e2_3,e_1,e_2を含む)
nde:=18: #要素自由度数

read "procs3.m":

# C;
# Cn;

Wall:=Vector(nde):
for i from 1 to nde do
  Wall[i]:= W[i]:
end do:
# ExportVector("out_C2.dat",C,target=delimited,delimiter=" "):
# ExportVector("out_Cn2.dat",Cn,target=delimited,delimiter=" "):
# ExportVector("out_W_elm2.dat",Wall,target=delimited,delimiter=" "):


# H:=diff 1次
H:=Array(1..nue,1..nde):
H2:=Array(1..nue,1..nde,1..nde):
H3:=Array(1..nue,1..nde,1..nde,1..nde):
H4:=Array(1..nue,1..nde,1..nde,1..nde,1..nde):

# Wzero:=Array[seq(1.0e-15,i=1..ndt)]:

# diff 1次 (Gamma matrix)
for i from 1 to nue do
  for j from 1 to nde do
    #2次以上の項がある場合はmtaylorが必要
    fff:=simplify(diff(C[i],W[j]));
    tmpeq:=normal(mtaylor(fff,W,tyl));
    save tmpeq, cat("./H/H_",i,"_",j,".m"):
    save tmpeq, cat("./H/H_",i,"_",j,".mpl"):
    H[i,j]:=tmpeq;
  end do;
end do;

save H, "./H/H.m":
save H, "./H/H.mpl":

# read "./H/H.m":

if ord>1 then
# diff 2次 （対称性を利用）
for i from 1 to nue do
  for j from 1 to nde do
    for k from j to nde do
      tmpeq:=diff(H[i,j],W[k]);
      save tmpeq, cat("./H2/H2_",i,"_",j,"_",k,".m"):
      save tmpeq, cat("./H2/H2_",i,"_",j,"_",k,".mpl"):
      H2[i,j,k]:=tmpeq;
      save tmpeq, cat("./H2/H2_",i,"_",k,"_",j,".m"):
      save tmpeq, cat("./H2/H2_",i,"_",k,"_",j,".mpl"):
      H2[i,k,j]:=tmpeq;
    end do;
  end do;
end do;
end if;

save H2, "./H2/H2.m":
save H2, "./H2/H2.mpl":

# quit:

# read "./H2/H2.m":

if ord>=3 then
# diff 3次 （対称性を利用）
for i from 1 to nue do
  for j from 1 to nde do
    for k from j to nde do
      for m from k to nde do
        tmpeq:=diff(H2[i,j,k],W[m]);
        save tmpeq, cat("./H3/H3_",i,"_",j,"_",k,"_",m,".m"):
        save tmpeq, cat("./H3/H3_",i,"_",j,"_",k,"_",m,".mpl"):
        H3[i,j,k,m]:=tmpeq;
        save tmpeq, cat("./H3/H3_",i,"_",j,"_",m,"_",k,".m"):
        save tmpeq, cat("./H3/H3_",i,"_",j,"_",m,"_",k,".mpl"):
        H3[i,j,m,k]:=tmpeq;
        save tmpeq, cat("./H3/H3_",i,"_",k,"_",j,"_",m,".m"):
        save tmpeq, cat("./H3/H3_",i,"_",k,"_",j,"_",m,".mpl"):
        H3[i,k,j,m]:=tmpeq;
        save tmpeq, cat("./H3/H3_",i,"_",k,"_",m,"_",j,".m"):
        save tmpeq, cat("./H3/H3_",i,"_",k,"_",m,"_",j,".mpl"):
        H3[i,k,m,j]:=tmpeq;
        save tmpeq, cat("./H3/H3_",i,"_",m,"_",j,"_",k,".m"):
        save tmpeq, cat("./H3/H3_",i,"_",m,"_",j,"_",k,".mpl"):
        H3[i,m,j,k]:=tmpeq;
        save tmpeq, cat("./H3/H3_",i,"_",m,"_",k,"_",j,".m"):
        save tmpeq, cat("./H3/H3_",i,"_",m,"_",k,"_",j,".mpl"):
        H3[i,m,k,j]:=tmpeq;
      end do;
    end do;
  end do;
end do;
end if;

save H3, "./H3/H3.m":
save H3, "./H3/H3.mpl":

# read "./H3/H3.m":

#ここから未修正171228
# if ord>=4 then
# # diff 4次（対称性を利用）
# for i from 1 to nu do
#   for j from 1 to ndt do
#     for k from j to ndt do
#       for m from k to ndt do
#         for s from m to ndt do
#           H4[i,j,k,m,s]:=diff(H3[i,j,k,m],W[s]);
#         end do;
#       end do;
#     end do;
#   end do;
# end do;
# end if;

##################################


# W0での評価
H:=subs(U1[1]=0.0,U1[2]=0.0,U1[3]=0.0,Tht1[1]=0.0,Tht1[2]=0.0,Tht1[3]=0.0,
U2[1]=0.0,U2[2]=0.0,U2[3]=0.0,Tht2[1]=0.0,Tht2[2]=0.0,Tht2[3]=0.0,
V[1]=0.0,V[2]=0.0,V[3]=0.0,Ps[1]=0.0,Ps[2]=0.0,Ps[3]=0.0,H):

H2:=subs(U1[1]=0.0,U1[2]=0.0,U1[3]=0.0,Tht1[1]=0.0,Tht1[2]=0.0,Tht1[3]=0.0,
U2[1]=0.0,U2[2]=0.0,U2[3]=0.0,Tht2[1]=0.0,Tht2[2]=0.0,Tht2[3]=0.0,
V[1]=0.0,V[2]=0.0,V[3]=0.0,Ps[1]=0.0,Ps[2]=0.0,Ps[3]=0.0,H2):

H3:=subs(U1[1]=0.0,U1[2]=0.0,U1[3]=0.0,Tht1[1]=0.0,Tht1[2]=0.0,Tht1[3]=0.0,
U2[1]=0.0,U2[2]=0.0,U2[3]=0.0,Tht2[1]=0.0,Tht2[2]=0.0,Tht2[3]=0.0,
V[1]=0.0,V[2]=0.0,V[3]=0.0,Ps[1]=0.0,Ps[2]=0.0,Ps[3]=0.0,H3):

save H, "./H/H_a.mpl":
save H2, "./H2/H2_a.mpl":
save H3, "./H3/H3_a.mpl":
save H, "./H/H_a.m":
save H2, "./H2/H2_a.m":
save H3, "./H3/H3_a.m":

quit;

# read "./H/H_a.m":
# read "./H2/H2_a.m":
# read "./H3/H3_a.m":



for i from 1 to nu do
  for j from 1 to ndt do
    H[i,j]:=evalf(H[i,j]);
    if ord>1 then
      for k from j to ndt do
        H2[i,j,k]:=evalf(H2[i,j,k]);
        H2[i,k,j]:=H2[i,j,k];
        if ord>=3 then
          for m from k to ndt do
            H3[i,j,k,m]:=evalf(H3[i,j,k,m]);
            H3[i,j,m,k]:=H3[i,j,k,m];
            H3[i,k,j,m]:=H3[i,j,k,m];
            H3[i,k,m,j]:=H3[i,j,k,m];
            H3[i,m,j,k]:=H3[i,j,k,m];
            H3[i,m,k,j]:=H3[i,j,k,m];
            for s from m to ndt do
              H4[i,j,k,m,s]:=evalf(H4[i,j,k,m,s]);
              H4[i,j,k,s,m]:=H4[i,j,k,m,s];
              H4[i,j,m,k,s]:=H4[i,j,k,m,s];
              H4[i,j,m,s,k]:=H4[i,j,k,m,s];
              H4[i,j,s,k,m]:=H4[i,j,k,m,s];
              H4[i,j,s,m,k]:=H4[i,j,k,m,s];
              H4[i,k,j,m,s]:=H4[i,j,k,m,s];
              H4[i,k,j,s,m]:=H4[i,j,k,m,s];
              H4[i,k,m,j,s]:=H4[i,j,k,m,s];
              H4[i,k,m,s,j]:=H4[i,j,k,m,s];
              H4[i,k,s,j,m]:=H4[i,j,k,m,s];
              H4[i,k,s,m,j]:=H4[i,j,k,m,s];
              H4[i,m,j,k,s]:=H4[i,j,k,m,s];
              H4[i,m,j,s,k]:=H4[i,j,k,m,s];
              H4[i,m,k,j,s]:=H4[i,j,k,m,s];
              H4[i,m,k,s,j]:=H4[i,j,k,m,s];
              H4[i,m,s,j,k]:=H4[i,j,k,m,s];
              H4[i,m,s,k,j]:=H4[i,j,k,m,s];
              H4[i,s,j,k,m]:=H4[i,j,k,m,s];
              H4[i,s,j,m,k]:=H4[i,j,k,m,s];
              H4[i,s,k,j,m]:=H4[i,j,k,m,s];
              H4[i,s,k,m,j]:=H4[i,j,k,m,s];
              H4[i,s,m,j,k]:=H4[i,j,k,m,s];
              H4[i,s,m,k,j]:=H4[i,j,k,m,s];
            end do;
          end do;
        end if;
      end do;
    end if;
end do;

tmp:=Matrix(nu,ndt):
HH:=Matrix(nu,ndf):
if ord>=2 then
  HH2:=Array(1..nu,1..ndf,1..ndf):
end if;
if ord>=3 then
  HH3:=Array(1..nu,1..ndf,1..ndf,1..ndf):
end if;

if ord>=4 then
  HH4:=Array(1..nu,1..ndf,1..ndf,1..ndf,1..ndf):
end if;

for i from 1 to nu do
  for j from 1 to ndf do
    HH[i,j]:=H[i,idx[j]]:
    tmp[i,j]:=H[i,j]:
    if ord>=2 then
      for k from 1 to ndf do
        HH2[i,j,k]:=Re(H2[i,idx[j],idx[k]]):
        if ord>=3 then
          for m from 1 to ndf do
            HH3[i,j,k,m]:=Re(H3[i,idx[j],idx[k],idx[m]]);
            if ord>=4 then
              for s from 1 to ndf do
                HH4[i,j,k,m,s]:=Re(H4[i,idx[j],idx[k],idx[m],idx[s]]);
              end do;
            end if;
          end do;
        end if;
      end do;
    end if;
  end do;
end do;

# ファイルへの出力
if ord>=2 then
  tmp2:=Matrix(ndf*ndf,nu):
end if;
if ord>=3 then
  tmp3:=Matrix(ndf*ndf*ndf,nu):
end if;

if ord>=4 then
  tmp4:=Matrix(ndf*ndf*ndf*ndf,nu):
end if;

for i from 1 to nu do
  for j from 1 to ndf do
    if ord>1 then
      for k from 1 to ndf do
        tmp2[(j-1)*ndf+k,i]:=HH2[i,j,k]:
        if ord>=3 then
          for m from 1 to ndf do
            tmp3[(j-1)*ndf*ndf+(k-1)*ndf+m,i]:=HH3[i,j,k,m]:
            if ord>=4 then
              for s from 1 to ndf do
                tmp4[(j-1)*ndf*ndf*ndf+(k-1)*ndf*ndf+(m-1)*ndf+s,i]:=HH4[i,j,k,m,s]:
              end do;
            end if;
          end do;
        end if;
      end do;
    end if;
  end do;
end do;

ExportMatrix("HH.dat",HH,delimiter=" ");
ExportMatrix("H.dat",tmp,delimiter=" ");
if ord>=2 then
  ExportMatrix("HH2.dat",tmp2,delimiter=" ");
end if;
if ord>=3 then
  ExportMatrix("HH3.dat",tmp3,delimiter=" ");
end if;
if ord>=4 then
  ExportMatrix("HH4.dat",tmp4,delimiter=" ");
end if;


# 特異値分解

r:=Rank(HH);
S, U, Vt:=SingularValues(HH,output=['S','U','Vt']):

ExportVector("out_sv.dat",S,delimiter=" ");
#ExportVector("out_left.dat",Column(U,nu),delimiter=" ");
#ExportVector("out_right.dat",Transpose(Row(Vt,ndf)),delimiter=" ");

p:=ndf - r;
q:=nu - r;

Hmat:=Matrix(ndf,p):
for i from 1 to p do
  for j from 1 to ndf do
    Hmat[j,i]:= Vt[ndf-(i-1),j]:
  end do:
end do:

ExportMatrix("out_Hmat.dat",Hmat,delimiter=" ");


print("---------RESULT-----------");

p;
q;
r;

quit;
