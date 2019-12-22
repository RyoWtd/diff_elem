# diff_elem
ヒンジ骨組の高次非適合成分計算プログラム

diff_main.mpl

インプットファイル説明

inp_param.csv
 モデル概要データ
 順に
 nu0 : 部材数
 nd0 : 節点数
 ndr : 固定自由度数
 nc : 追加拘束条件数

inp_fix.csv
 固定自由度定義データ
 各行 1自由度ずつ定義する
 1列目 : 節点番号
 2列目 : 固定自由度（1～6）

inp_coord.csv
 節点座標データ
 各行 節点番号順に X,Y,Z座標を指定（数式入力可）

inp_cnn.csv
 部材接続関係定義データ
 部材数(nu0)と同数行を部材番号順に定義
 1列目 : 部材始端節点番号
 2列目 : 部材終端節点番号

inp_n_hinge.csv
 部材端ヒンジ数定義データ
 部材数(nu0)と同数行を部材番号順に定義
 1列目 : 始端ヒンジ軸数（1:ヒンジ, 2:ユニバーサルジョイント）
 2列目 : 終端ヒンジ軸数（1:ヒンジ, 2:ユニバーサルジョイント）

inp_hinge_ind.csv
 追加ヒンジの軸ベクトル設定位置定義データ
 追加ヒンジ軸数と同行数を定義する
 1列目 : 定義要素番号
 2列目 : 定義位置 始終端指定(1 or 2)
 3列目 : 同位置での軸番号指定
 (1軸のみなら1のみ, ユニバーサルジョイントなら1,2それぞれ1行ずつ指定する)

inp_hinge_vec.csv
 追加ヒンジの軸ベクトル成分定義データ
 追加ヒンジ軸と同行数を定義する（順番はinp_hinge_ind.datに対応）
 1,2,3列目 : ヒンジ軸ベクトル成分（数式入力可）
 (要素座標系で定義する。大きさ1で正規化しておく)

