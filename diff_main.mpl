with(LinearAlgebra):

# 計算次数
ord:=3:
tyl:=3:

# inp_param.dat
#  モデル概要データ
#  順に
#  nu0 : 部材数
#  nd0 : 節点数
#  ndr : 固定自由度数
#  nc : 追加拘束条件数
inpparam:=ImportVector("inp_param.csv"):
nu0:=inpparam[1]:
nd0:=inpparam[2]:
ndr:=inpparam[3]:
nc:=inpparam[4]:
# 全自由度数
ndt:=nd0*6+nu0*6:

nue:=20: #要素非適合量数
# (12成分+ 8成分 e1_1,e2_1,e3_1,e2_1,e2_2,e2_3,e_1,e_2を含む)
# この中から選ぶ
nde:=18: #要素自由度数

nd:=ndt-ndr: #固定自由度数を除いた全自由度数
# 全体自由度の順番:はじめに全節点の並進UU*3,回転TT*3,
#                次に各要素中央点の並進VV*3,回転PP*3

# 座標定義
# inp_coord.csvから読み込み
xvg:=[seq(Vector([0,0,0]),i=1..nd0)]:
temp:=ImportMatrix("inp_coord.csv", delimiter=","):
for i from 1 to nd0 do
  for j from 1 to 3 do
    if type(temp[i][j],string) = true then
      xvg[i][j]:=parse(temp[i][j]):
    else
      xvg[i][j]:=temp[i][j]:
    end if:
  end do:
end do:
# xvg[2][1]:=1.0:
# xvg[3][1]:=2.0:

# idx: 固定する自由度
# 固定自由度をinp_fix.csvから読み込み
# idxf:=Vector([1,2,3,4,6,14,15]): #問題設定
idx:=Vector(ndt):
for i from 1 to ndt do
    idx[i]:=i:
end do:
if ndr > 0 then
    temp2:=ImportMatrix("inp_fix.csv", delimiter=","):
    idxf:=Vector(ndr):
    for i from 1 to ndr do
        idxf[i]:= (temp2[i][1]-1)*6 + temp2[i][2]:
    end do:
    for i from 1 to ndr do
        idx[idxf[i]] := 0:
    end do:
    k:=1:
    for i from 1 to ndt do
        if idx[i] <> 0 then
            idx[k] := idx[i]:
            k:=k+1:
        end if:
    end do:
    for i from k to ndt do
        idx[i] := 0:
    end do:
end if:

# 接続情報
#  (inp_cnn.dat から読み込み)
#  cnn    : 各部材の始終端節点番号
cnn:=[seq([0,0],i=1..nu0)]:
tempM:=ImportMatrix("inp_cnn.csv",delimiter=","):
for i from 1 to nu0 do
 for j from 1 to 2 do
  cnn[i][j] := tempM[i,j]:
 end do:
end do:
# cnn[1][1]:=1:
# cnn[1][2]:=2:
# cnn[2][1]:=2:
# cnn[2][2]:=3:

# 回転ヒンジの数情報
# 1:1本（レボリュートジョイント）, 2:2本（ユニバーサルジョイント）
#  (inp_n_hinge.dat から読み込み)
#  ih     : 各部材の始終端ヒンジ軸数（1:ヒンジ,2:ユニバーサルジョイント)
ih:=[seq([0,0],i=1..nu0)]:
tempM:=ImportMatrix("inp_n_hinge.csv",delimiter=","):
for i from 1 to nu0 do
 for j from 1 to 2 do
  ih[i][j] := tempM[i,j]:
 end do:
end do:
# ih[2][1]:=1:

# 局所座標でのヒンジ方向
# 1箇所に2軸まで定義可能とした

# ヒンジ方向の読み込み（各部材局所座標系で指定）
#  (inp_hinge_ind.dat, inp_hinge_vec.dat より読み込み）
#  hht[i][j][k]  : 部材iのj端側k番目のヒンジ軸の方向ベクトル（局所座標系）
#  hhg[i][j][k] : 部材iのj端側k番目のヒンジ軸の方向ベクトル（全体座標系）
hht:=[seq[seq(seq(Vector([0,0,0]),l=1..2),j=1..2)],i=1..nu0]:

# hht[2][1][1]:=Vector([sqrt(2)/2,sqrt(2)/2,0]): # model B
# hht[2][1][1]:=Vector([0,1,0]): # model A

tempM2:=ImportMatrix("inp_hinge_ind.csv", delimiter=","):
tempM3:=ImportMatrix("inp_hinge_vec.csv", delimiter=","):
aa:=RowDimension(tempM2);
tempM3temp:=Vector(3):
for i from 1 to aa do
  for j from 1 to 3 do
    if type(tempM3[i,j],string) = true then
      tempM3temp[j] := parse(tempM3[i,j]):
    else
      tempM3temp[j] := tempM3[i,j]:
    end if:
  end do:
  hht[tempM2[i,1]][tempM2[i,2]][tempM2[i,3]]:=Vector([tempM3temp[1],tempM3temp[2],tempM3temp[3]]):
end do:




# IDM[k] : 削減前自由度番号kの削減後自由度番号
IDM:=Vector(ndt):
ii:=1:
jj:=0:
if ndr > 0 then
    for i from 1 to ndt do
        if ii<=ndr then
            if i=idxf[ii] then
                IDM[i]:=0:
                ii:=ii+1:
            else
                jj:=jj+1:
                IDM[i]:=jj:
            end if:
        else
            jj:=jj+1:
            IDM[i]:=jj:
        end if;
    end do;
else
    for i from 1 to ndt do
        jj:=jj+1:
        IDM[i]:=jj:
    end do:
end if:
ExportVector("IDM.dat",IDM,target=delimited,delimiter=" "):

# 要素自由度番号-全体自由度番号（固定自由度削減前）関係の作成
LMOF:=[seq(Vector(nde),j=1..nu0)]:
# 全体自由度の順番:はじめに全節点の並進UU*3,回転TT*3,
#                次に各要素中央点の並進VV*3,回転PP*3
# 要素自由度の順番:はじめに始端点の並進UU*3,回転TT*3,
#                次に終端点の並進UU*3,回転TT*3,
#                最後に要素中央点の並進VV*3,回転PP*3
# LMOF[j][k] : 要素jのk番目自由度の全体自由度中の番号
for j from 1 to nu0 do
    kkk1:=cnn[j][1]:
    kkk2:=cnn[j][2]:
    for k from 1 to 6 do
        LMOF[j][k]:=6*(kkk1-1)+k:
        LMOF[j][k+6]:=6*(kkk2-1)+k:
    end do;
    for k from 1 to 6 do
        LMOF[j][k+12]:=6*nd0+6*(j-1)+k:
    end do;
end do;

# 要素非適合成分番号-全体非適合成分番号関係の作成
LMOC:=[seq(Vector(nue),j=1..nu0)]:
# 全体非適合成分番号の順番:各要素の非適合成分番号を要素順に並べる
# 要素非適合成分番号の順番:
#   DeltaU_11~13,DeltaU_21~23,DeltaTht_11~13,DeltaTht_21~23,
#   e1_1~3,e2_1~3,e_1,e_2 の順
# LMOC[j][k] : 要素jのk番目非適合成分番号の全体自由度中の番号
#              採用しない場合は0を格納する

ii:=0:
for j from 1 to nu0 do
    # DeltaU_11~13,DeltaU_21~23
    for k from 1 to 6 do
        ii:=ii+1:
        LMOC[j][k]:=ii:
    end do;
    # DeltaTht_11~13
    for k from 1 to 3 do
        if ih[j][1]=0 then
            ii:=ii+1:
            LMOC[j][k+6]:=ii:
        else
            LMOC[j][k+6]:=0:
        end if:
    end do;
    # DeltaTht_21~23
    for k from 1 to 3 do
        if ih[j][2]=0 then
            ii:=ii+1:
            LMOC[j][k+9]:=ii:
        else
            LMOC[j][k+9]:=0:
        end if:
    end do;
    # e1_1~3
    if ih[j][1]=1 then
        inswitch := 3: #仮対処　後日正しく選択するようになおす
        if inswitch = 3 then
            ii:=ii+1:
            LMOC[j][13]:=ii:
            ii:=ii+1:
            LMOC[j][14]:=ii:
            LMOC[j][15]:=0:
        elif inswitch = 2 then
            ii:=ii+1:
            LMOC[j][13]:=ii:
            LMOC[j][14]:=0:
            ii:=ii+1:
            LMOC[j][15]:=ii:
        elif inswitch = 1 then
            LMOC[j][13]:=0:
            ii:=ii+1:
            LMOC[j][14]:=ii:
            ii:=ii+1:
            LMOC[j][15]:=ii:
        end if:
    else
        for k from 1 to 3 do
                LMOC[j][k+12]:=0:
        end do;
    end if:
    # e2_1~3
    if ih[j][2]=1 then
        inswitch := 3: #仮対処　後日正しく選択するようになおす
        if inswitch = 3 then
            ii:=ii+1:
            LMOC[j][16]:=ii:
            ii:=ii+1:
            LMOC[j][17]:=ii:
            LMOC[j][18]:=0:
        elif inswitch = 2 then
            ii:=ii+1:
            LMOC[j][16]:=ii:
            LMOC[j][17]:=0:
            ii:=ii+1:
            LMOC[j][18]:=ii:
        elif inswitch = 1 then
            LMOC[j][16]:=0:
            ii:=ii+1:
            LMOC[j][17]:=ii:
            ii:=ii+1:
            LMOC[j][18]:=ii:
        end if:
    else
        for k from 1 to 3 do
                LMOC[j][k+15]:=0:
        end do;
    end if:
    # e1
    if ih[j][1]=2 then
        ii:=ii+1:
        LMOC[j][19]:=ii:
    else
        LMOC[j][19]:=0:
    end if:
    # e2
    if ih[j][2]=2 then
        ii:=ii+1:
        LMOC[j][20]:=ii:
    else
        LMOC[j][20]:=0:
    end if:
end do;
# 全非適合成分数
nu:=ii:

# 要素H行列を読み込み
read "./H/H_a.m":
read "./H2/H2_a.m":
read "./H3/H3_a.m":

# 要素行列を全体行列に組み込み
HH:=Matrix(nu,nd):
for i from 1 to nu0 do
    tmpH:=subs(
        XX1[1]=xvg[cnn[i][1]][1],
        XX1[2]=xvg[cnn[i][1]][2],
        XX1[3]=xvg[cnn[i][1]][3],
        XX2[1]=xvg[cnn[i][2]][1],
        XX2[2]=xvg[cnn[i][2]][2],
        XX2[3]=xvg[cnn[i][2]][3],
        hvl11=hht[i][1][1],
        hvl12=hht[i][1][2],
        hvl21=hht[i][2][1],
        hvl22=hht[i][2][2],
        H):
    for j from 1 to nue do
        if LMOC[i][j] <> 0 then
            irow := LMOC[i][j]:
            for k from 1 to nde do
                icol := IDM[LMOF[i][k]]:#削減後の自由度番号
                if icol<>0 then
                    HH[irow,icol]:=HH[irow,icol]+rationalize(tmpH[j,k]):
                end if:
            end do;
        end if:
    end do;
end do;
HH:=simplify(HH):
save HH, "./out/HH.mpl";
save HH, "./out/HH.m";

HH2:=Array(1..nu,1..nd,1..nd):
for i from 1 to nu0 do
    tmpH2:=subs(
        XX1[1]=xvg[cnn[i][1]][1],
        XX1[2]=xvg[cnn[i][1]][2],
        XX1[3]=xvg[cnn[i][1]][3],
        XX2[1]=xvg[cnn[i][2]][1],
        XX2[2]=xvg[cnn[i][2]][2],
        XX2[3]=xvg[cnn[i][2]][3],
        hvl11=hht[i][1][1],
        hvl12=hht[i][1][2],
        hvl21=hht[i][2][1],
        hvl22=hht[i][2][2],
        H2):
    for j from 1 to nue do
        if LMOC[i][j] <> 0 then
            irow := LMOC[i][j]:
            for k from 1 to nde do
                icol1 := IDM[LMOF[i][k]]:#削減後の自由度番号
                if icol1<>0 then
                    for l from 1 to nde do
                        icol2 := IDM[LMOF[i][l]]:
                        if icol2<>0 then
                            HH2[irow,icol1,icol2]:=HH2[irow,icol1,icol2]
                            +rationalize(tmpH2[j,k,l]):
                        end if:
                    end do;
                end if:
            end do;
        end if:
    end do;
end do;
HH2:=simplify(HH2):
save HH2, "./out/HH2.mpl";
save HH2, "./out/HH2.m";

HH3:=Array(1..nu,1..nd,1..nd,1..nd):
for i from 1 to nu0 do
    tmpH3:=subs(
        XX1[1]=xvg[cnn[i][1]][1],
        XX1[2]=xvg[cnn[i][1]][2],
        XX1[3]=xvg[cnn[i][1]][3],
        XX2[1]=xvg[cnn[i][2]][1],
        XX2[2]=xvg[cnn[i][2]][2],
        XX2[3]=xvg[cnn[i][2]][3],
        hvl11=hht[i][1][1],
        hvl12=hht[i][1][2],
        hvl21=hht[i][2][1],
        hvl22=hht[i][2][2],
        H3):
    for j from 1 to nue do
        if LMOC[i][j] <> 0 then
            irow := LMOC[i][j]:
            for k from 1 to nde do
                icol1 := IDM[LMOF[i][k]]:#削減後の自由度番号
                if icol1<>0 then
                    for l from 1 to nde do
                        icol2 := IDM[LMOF[i][l]]:
                        if icol2<>0 then
                            for ll from 1 to nde do
                                icol3 := IDM[LMOF[i][ll]]:
                                if icol3<>0 then
                                    HH3[irow,icol1,icol2,icol3]:=
                                    HH3[irow,icol1,icol2,icol3]
                                    +rationalize(tmpH3[j,k,l,ll]):
                                end if:
                            end do;
                        end if:
                    end do;
                end if:
            end do;
        end if:
    end do;
end do;
HH3:=simplify(HH3):
save HH3, "./out/HH3.mpl";
save HH3, "./out/HH3.m";

# Octave用出力
param:=Vector(3):
param[1]:=nu:
param[2]:=nd:
param[3]:=nd+ndr:
ExportVector("./out/out_param.dat",param,delimiter=" "):

tmp2:=Matrix(nd*nd,nu):
tmp3:=Matrix(nd*nd*nd,nu):
for i from 1 to nu do
  for j from 1 to nd do
    if ord>1 then
      for k from 1 to nd do
        tmp2[(j-1)*nd+k,i]:=HH2[i,j,k]:
        if ord=3 then
          for m from 1 to nd do
            tmp3[(j-1)*nd*nd+(k-1)*nd+m,i]:=HH3[i,j,k,m]:
          end do;
        end if;
      end do;
    end if;
  end do;
end do;

ExportMatrix("./out/HH.dat",HH,delimiter=" ");
ExportMatrix("./out/HH2.dat",tmp2,delimiter=" ");
ExportMatrix("./out/HH3.dat",tmp3,delimiter=" ");

# 特異値分解

rrr:=Rank(HH);
S, U, Vt:=SingularValues(HH,output=['S','U','Vt']):

ExportVector("sv.dat",S,delimiter=" ");
ExportVector("left.dat",Row(U,nu),delimiter=" ");

ExportVector("right.dat",Transpose(Row(Vt,nd)),delimiter=" ");


quit;



# # read C, Cn, W
# read "procs.m"