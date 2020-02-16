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