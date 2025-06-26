#完事后在inferCNV的output里找到 HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.cell_groupings

#检查生成
ls "/Users/leishengyue/Desktop/data/论文数据/inferCNV/MPS1iclone4forUphy" | grep cell_groupings
#看见17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings

#进入文件夹
cd "/Users/leishengyue/Desktop/data/论文数据/inferCNV/MPS1iclone4forUphy"

# 删除引用/对照细胞行（以 all_references 开头），生成新的 trimmed 文件
sed '/^CTRL/d' "17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings" > MPS1iClone1.cell_groupings

#找一下Uphyplot2文件夹在哪
#在/Users/leishengyue/Desktop/data/UPhyloplot2 里

#确保input文件夹存在
mkdir -p /Users/leishengyue/Desktop/data/UPhyloplot2

#把cell_groupings复制进input里去
cp "/Users/leishengyue/Desktop/data/论文数据/inferCNV/MPS1iclone1forUphy/MPS1iClone1.cell_groupings" /Users/leishengyue/Desktop/data/UPhyloplot2/Inputs/

#运行
cd /Users/leishengyue/Desktop/data/UPhyloplot2
python3 uphyloplot2.py
