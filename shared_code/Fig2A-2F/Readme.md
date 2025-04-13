# Figure of Similarity bar (SIMILE.pdf)
原始是fig2， 改成 correspondence 后改成了 fig 1B-1G
使用吕华伟给的数据集验证效果

## 0.运行 GNPS library Search 和 MNA 
    Library 具体参数见文件夹下`params.xml`
    MNA 使用 modified cosine

## 1.Check if how much part these standards in EDB or ISDB
    `Check_in_db.py`

## 2.Filter the annotation results of GNPS and MNA
有些GNPS中的数据是没有smiles的
**手动补充没有 smiles 的 CCMS， 以及result中所有scan的smiles！！！！** 
Sync to `ccms补充信息的.xlsx`
可以直接查看 0.1 阈值中,CCMS缺失的smiles

### 2.1 `MatRes_fbmn.py`, `MatRes_mna.py`
    GNPS library search的结果不多的话，可以全部算出mcs再划分；
    设置 similarity, Peakpercentage , matched peaks 的 threshold筛选MNA产生的edbMS1, npMS1的结果文件
    EDB和ISDB分开处理，生成过滤后 match 的结果
    生成之后再确认一下有没有缺失smiles
    
    GNPS_LibSearch
    补充进对应 #Scan# 的smiles
    
### 2.2 `MCS_sta_fbmn.py`,`MCS_sta_mna.py`, `MCS_sta_mix.py`
    根据 SIMILE的 Max Common Structure算法和相似度标准，将 annotation 分为 5类
    High Similarity (>70%)
    Medium Similarity (45%-70%)
    Low Similarity (35%-45%)
    Not Similar (<35%)
    Not compared

### 2.3 Confusion matrix
            Not Selected (< Threshold)  Selected (>Threshold)
MCS < 0.15(0.35)           True Neg              False Pos
MCS > 0.7            False Neg             True Pos
`GNPS_confustion_matrix.py`和`MNA_confusion_matrix.py`
在GNPS中进行library Search，阈值： pair similarity >= 0.1, shared peaks >= 1, top 100
根据指定阈值将搜库结果分为df_selected和df_unselected,
scans in df_selected + scans in df_unselected = 99，要除去unselected中selected的ID
根据MCS的值得到对应的混淆矩阵
MCS < 0.15 和 MCS < 0.35 的算出来的结果会有些许差异

MNA的结果则更为复杂，依据一个规则：同一个feature的不同方法及不同数据库匹配的结果都取最大值
scans in df_selected + scans in df_unselected = 99，要除去unselected中selected的ID

### 4. 根据
