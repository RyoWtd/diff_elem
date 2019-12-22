# 更新履歴
# 2019.12.1
#   作成開始
#
#
#
# W[1]~W[3]:U1[1],U1[2],U1[3]
# W[4]~W[6]:Tht1[1],Tht1[2],Tht1[3]
# W[7]~W[9]:U2[1],U2[2],U2[3]
# W[10]~W[12]:Tht2[1],Tht2[2],Tht2[3]
# W[13]~W[15]:V[1],V[2],V[3]
# W[16]~W[18]:Ps[1],Ps[2],Ps[3]
#
# C[1]~C[3]:DeltaU_11, DeltaU_12, DeltaU_13
# C[4]~C[6]:DeltaU_21, DeltaU_22, DeltaU_23
# C[7]~C[9]:DeltaTht_11, DeltaTht_12, DeltaTht_13
# C[10]~C[12]:DeltaTht_21, DeltaTcdht_22, DeltaTht_23
# C[13]~C[15]:e1_1, e2_1, e3_1
# C[16]~C[18]:e1_2, e2_2, e3_2
# C[19]:e_1
# C[20]:e_2
#
# hvl11 : 始端でのl本目のヒンジ軸ベクトル(要素局所座標系)
# hvl12 : 始端での2本目のヒンジ軸ベクトル(要素局所座標系)
# hvl21 : 終端でのl本目のヒンジ軸ベクトル(要素局所座標系)
# hvl22 : 終端での2本目のヒンジ軸ベクトル(要素局所座標系)
# 要素局所座標は x:材軸方向, y:g-Z X l-x, z: x X y で定める
#
# xv[1]:始点座標
# xv[2]:終点座標
#
# 要素局所座標

# 計算に使用する仮数部桁数を変更
Digits:=18;

with(LinearAlgebra):

# ord:=3: #次数
# tyl:=3: #テイラー展開次数

nue:=20: #要素非適合量数
# (12成分+ 8成分 e1_1,e2_1,e3_1,e2_1,e2_2,e2_3,e_1,e_2を含む)
nde:=18: #要素自由度数

# 始端節点併進自由度
UU1:=Vector(3,symbol=U1):
# 始端節点回転自由度
TT1:=Vector(3,symbol=Tht1):
# 終端節点併進自由度
UU2:=Vector(3,symbol=U2):
# 終端節点回転自由度
TT2:=Vector(3,symbol=Tht2):
# 部材中央点併進自由度
VV:=Vector(3,symbol=V):
# 部材中央点回転自由度
PP:=Vector(3,symbol=Ps):

# 節点座標
xv[1]:=Vector(3,symbol=XX1): #xv[1]:始点
xv[2]:=Vector(3,symbol=XX2): #xv[2]:終点

#W:=Array[seq(0,i=1..nde)]:
for j from 1 to 3 do
  W[j]:=UU1[j]:
  W[3+j]:=TT1[j]:
  W[6+j]:=UU2[j]:
  W[9+j]:=TT2[j]:
  W[12+j]:=VV[j]:
  W[15+j]:=PP[j]:
end do:

Wall:=Vector(nde):
for i from 1 to nde do
  Wall[i]:= W[i]:
end do:
ExportVector("out_W_elm.dat",Wall,target=delimited,delimiter=" "):

# 部材長
LEN:=sqrt((xv[2][1]-xv[1][1])**2 + (xv[2][2]-xv[1][2])**2
      +(xv[2][3]-xv[1][3])**2 ):

zz:=Vector([0,0,1]):
t1:=Vector(3):
t2:=Vector(3):
t3:=Vector(3):

# 部材局所座標
tv:=[seq(Vector(3),i=1..3)]:
tv[1]:=(xv[2]-xv[1])/LEN:
t1:=tv[1]:
t2:=CrossProduct(zz,t1):
t3:=CrossProduct(t1,t2):
tv[2]:=t2:
tv[3]:=t3:

# 局所座標でのヒンジ方向
# 1箇所に二軸まで追加できるように修正
hh11:=Vector(3,symbol=hvl11):
hh12:=Vector(3,symbol=hvl12):
hh21:=Vector(3,symbol=hvl21):
hh22:=Vector(3,symbol=hvl22):
# hvl11 : 始端でのl本目のヒンジ軸ベクトル
# hvl12 : 始端での2本目のヒンジ軸ベクトル
# hvl21 : 終端でのl本目のヒンジ軸ベクトル
# hvl22 : 終端での2本目のヒンジ軸ベクトル

# 全体座標でのヒンジ方向
hhg:=[seq(Vector(3),l=1..2)]:
hhg[1][1]:=hh11[1]*tv[1]+hh11[2]*tv[2]+hh11[3]*tv[3];
hhg[1][2]:=hh12[1]*tv[1]+hh12[2]*tv[2]+hh12[3]*tv[3];
hhg[2][1]:=hh21[1]*tv[1]+hh21[2]*tv[2]+hh21[3]*tv[3];
hhg[2][2]:=hh22[1]*tv[1]+hh22[2]*tv[2]+hh22[3]*tv[3];


# 部材回転角
# nn:=Vector(3):
phi:=VectorNorm(PP,2,conjugate=false):
nn:=PP/phi:

nnn:=Vector(3):

# 変形後のヒンジ方向（部材側）
hhg11:=[seq(Vector(3),i=1..2)]:
hhg12:=[seq(Vector(3),i=1..2)]:
for j from 1 to 2 do
  # 1番目の軸ベクトルを回転する
  tt:=hhg[j][1]:
  nnn:=nn*DotProduct(nn,tt,conjugate=false):
  # 三角関数 5次項まで
  # hhg11[j]:=nnn + (1-phi**2/2+phi**4/24)*(tt-nnn)
  #   - CrossProduct(tt,nn)*(phi-phi**3/6+phi**5/120):
  # 三角関数 4次項まで
  # hhg11[j]:=nnn + (1-phi**2/2+phi**4/24)*(tt-nnn)
  #   - CrossProduct(tt,nn)*(phi-phi**3/6):
  # 三角関数 3次項まで
  hhg11[j]:=nnn + (1-phi**2/2)*(tt-nnn)
    - CrossProduct(tt,nn)*(phi-phi**3/6):
  # 2番目の軸ベクトルを回転する
  tt:=hhg[j][2]:
  nnn:=nn*DotProduct(nn,tt,conjugate=false):
  # 三角関数 5次項まで
  # hhg12[j]:=nnn + (1-phi**2/2+phi**4/24)*(tt-nnn)
  #   - CrossProduct(tt,nn)*(phi-phi**3/6+phi**5/120):
  # 三角関数 4次項まで
  # hhg12[j]:=nnn + (1-phi**2/2+phi**4/24)*(tt-nnn)
  #   - CrossProduct(tt,nn)*(phi-phi**3/6):
  # 三角関数 3次項まで
  hhg12[j]:=nnn + (1-phi**2/2)*(tt-nnn)
    - CrossProduct(tt,nn)*(phi-phi**3/6):
end do:

# 節点回転角
phi2:=Vector(2):
phi2[1]:=VectorNorm(TT1,2,conjugate=false):
phi2[2]:=VectorNorm(TT2,2,conjugate=false):
nn2[1]:=TT1/phi2[1]:
nn2[2]:=TT2/phi2[2]:

# 変形後のヒンジ方向（節点側）
hhg2:=[seq(Vector(3),i=1..2)]:
for j from 1 to 2 do #1:始端側 2:終端側
  tt:=hhg[j][1];
  # いずれも1番目の軸ベクトルを回転する
  nnn:=nn2[j]*DotProduct(nn2[j],tt,conjugate=false);
  # 三角関数 5次項まで
  # hhg2[j]:=nnn + (1-phi2[j]**2/2+phi2[j]**4/24)*(tt-nnn)
  #   - CrossProduct(tt,nn2[j])*(phi2[j]-phi2[j]**3/6+phi2[j]**5/120);
  # 三角関数 4次項まで
  # hhg2[j]:=nnn + (1-phi2[j]**2/2+phi2[j]**4/24)*(tt-nnn)
  #   - CrossProduct(tt,nn2[j])*(phi2[j]-phi2[j]**3/6);
  # 三角関数 3次項まで
  hhg2[j]:=nnn + (1-phi2[j]**2/2)*(tt-nnn)
    - CrossProduct(tt,nn2[j])*(phi2[j]-phi2[j]**3/6);
end do:

# 変形後の部材局所座標
ttv:=[seq(Vector(3),i=1..3)]:
for j from 1 to 3 do
  tt:=tv[j]:
  nnn:=nn*DotProduct(nn,tt,conjugate=false):
  # 三角関数 5次項まで
  # ttv[j]:=nnn + (1-phi**2/2+phi**4/24)*(tt-nnn)
  #   - CrossProduct(tt,nn)*(phi-phi**3/6+phi**5/120):
  # 三角関数 4次項まで
  # ttv[j]:=nnn + (1-phi**2/2+phi**4/24)*(tt-nnn)
  #   - CrossProduct(tt,nn)*(phi-phi**3/6):
  # 三角関数 3次項まで
  ttv[j]:=nnn + (1-phi**2/2)*(tt-nnn)
    - CrossProduct(tt,nn)*(phi-phi**3/6):
end do:

# 材端変位の不適合量（並進・回転成分（ヒンジ無し端部））
# 全体座標表示
C1:=[seq(Vector(3),j=1..2)]:
C2:=[seq(Vector(3),j=1..2)]:
for j from 1 to 2 do
  if j=1 then
    C1[j]:=UU1 -VV -1/2*LEN*(tv[1]-ttv[1]);
    C2[j]:=TT1-PP;
  else
    C1[j]:=UU2 -VV +1/2*LEN*(tv[1]-ttv[1]);
    C2[j]:=TT2-PP;
  end if;
end do;

# 不適合量のベクトル
x1:=Vector(3):
x2:=Vector(3):
fff:=Vector(3):
C:=Vector(nue):
Cn:=Vector(nue):

ii:=0:


for j from 1 to 2 do
# DeltaU_jk 作成
  for k from 1 to 3 do
    ii:=ii+1:
    C[ii]:=C1[j][k]:
    Cn[ii]:= cat("DeltaU_",j,k):
  end do:
end do:
# DeltaTht_jk 作成(ヒンジ無し)
for j from 1 to 2 do
  for k from 1 to 3 do
    ii:=ii+1:
    C[ii]:=C2[j][k]:
    Cn[ii]:= cat("DeltaTht_",j,k):
  end do:
end do:
# e_kj 作成(ヒンジ1本)
for j from 1 to 2 do
  ddd:=CrossProduct(hhg11[j],hhg2[j]):
  for k from 1 to 3 do
    tt[k]:= ddd[k]
  end do:
  C[ii+1]:=simplify(expand(tt[1])):
  C[ii+2]:=simplify(expand(tt[2])):
  C[ii+3]:=simplify(expand(tt[3])):
  Cn[ii+1]:=cat("e1_",j):
  Cn[ii+2]:=cat("e2_",j):
  Cn[ii+3]:=cat("e3_",j):
  ii:=ii+3:
end do:
# e_j 作成(ヒンジ2本)
for j from 1 to 2 do
  C[ii+1]:=DotProduct(hhg12[j],hhg2[j],conjugate=false):
  Cn[ii+1]:=cat("e_",j):
  ii:=ii+1:
end do:

templl:=Vector(3):

ExportVector("out_C.dat",C,target=delimited,delimiter=" "):
ExportVector("out_Cn.dat",Cn,target=delimited,delimiter=" "):

save C, Cn, W, "procs3.m":

quit:
