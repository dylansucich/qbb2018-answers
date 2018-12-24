#!/Users/cmdb/miniconda3/bin/python 	

import scanpy.api as sc
import sys
import matplotlib
sc.settings.autoshow = False




adata = sc.read_10x_h5(sys.argv[1])
adata.var_names_make_unique()

sc.tl.pca(adata, svd_solver='auto')
sc.pl.pca(adata, save="PCA_adata.png", title="PCA Before Zheng Pre-processing")
# sc.tl.pca(adata, n_comps=50)

# var_names1 = adata.var_names
# for i in var_names1:
# 	print(adata.var_names[i])
	# if adata.var_names[str(var_names1[i])]:
	# 	print(adata.var_names[i])

# sc.pl.pca(adata, save="PCA_adata.png")

# sc.tl.tsne(adata)
# sc.pl.tsne(adata, save="tSNE_adata.png")

# sc.pp.neighbors(adata)
# sc.tl.umap(adata)
# sc.pl.umap(adata, save="uMAP_adata.png")


# filtering_method = ["recipe_zheng17", "recipe_weinreb17", "recipe_seurat"]
# var_names1 = adata.var_names["Slc17A7", "Slc17a6", "Gad1", "Gad2", "Slc32a1", "Fox3", "Gfap", "Sst", "Htr3a"]
# var_names1 = adata.var_names[0:]
nvars = ['Snap25', 'Syt1', 'Slc17A6', 'Slc32a1', 'Olig1', 'Sox9', 'Cldn5', 'C1qa']

sc.pp.recipe_zheng17(adata, n_top_genes=1000, log=True, plot=False, copy=False)

sc.tl.pca(adata, svd_solver='auto')
sc.pl.pca(adata, save="PCA_adata_Zheng.png", title="Zheng filtered PCA")

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=None)
sc.tl.louvain(adata)

sc.tl.umap(adata)
sc.pl.umap(adata, color="louvain", save="uMAP_adata_Zheng.png", title="uMAP of Zheng filtered scRNA-seq", legend_loc="right margin")

sc.tl.tsne(adata)
sc.pl.tsne(adata, color="louvain", save="tSNE_adata_Zheng.png", title="t-SNE of Zheng filtered scRNA-seq", legend_loc="right margin")

sc.tl.rank_genes_groups(adata, groupby="louvain", method="logreg")
sc.pl.rank_genes_groups(adata, save="Rank_GeneGroups_adata_Zheng_logreg.png")

sc.tl.rank_genes_groups(adata, groupby="louvain", method="t-test")
sc.pl.rank_genes_groups(adata, save="Rank_GeneGroups_adata_Zheng_t-test.png")


genes = adata.var.gene_ids.index.get_values()
GOI = ['Sst', 'Pdyn', 'Opkr1', 'Slc17a7','Slc17a6','Slc17a8', 'Slc32a1','Gad1','Gad2', 'Crhbp','Npy', 'Tac1','Tac2','Tph1','Tph2','Gfap','Pvalb', 'Syt2', 'Htr3a', 'Drd1', 'Drd2', 'Crh','Slc6a4','Slc6a3','Slc6a2','Cck','Oxt','Prkcd', 'Ddc','Pnmt','Pomc','Kcnmb1','Kcnj3','Kcnj6','Kcnj9','Kcnj5','Rbfox3','Fos','Jun', 'Lhx6','Cacna2d2','Erbb4','4932438A13Rik','Snhg11','Plec','Cdc42ep2','Htt','App', 'Ranbp2','Tox3','Nrip3','Slc7a3','Npy2r','Nxph2','Nos1','Camk2n2','Adcy2','Kcnip1','Ank1','Sox6', 'Sox17', 'Chat','Arx', 'Pax6'] 

present = [x for x in GOI if x in genes]
print(present)

present = ['Sst', 'Gad1', 'Gad2', 'Npy', 'Tac1', 'Tac2', 'Htr3a', 'Slc6a4', 'Cck', 'Fos', 'Jun', 'Lhx6', 'Cacna2d2', 'Erbb4', 'Snhg11', 'Nxph2', 'Sox17', 'Arx', 'Pax6', 'Olig1']

# 1 - Gad1 Gad2 Lhx6
# 5 - Tac2 Cck
# 14 - HTR3A Erbb4
# 17 - 
# 19 - Sox17


cluster_names = ['0', 'Sst GABAergic Interneurons', '2', 'Glial cells', '4', 'Tac2 Cck GABAergic Interneurons', '6', 'Radial glia 1', '8', '9', 'Radial glia 2', '11', '12', '13', 'Vip, Htr3a, Erbb4 GABAergic interneurons', '15', '16', 'GABAergic Interneurons', '18', 'Sox17 Oligodendrocytes', '20', '21', '22']

adata.rename_categories('louvain', cluster_names)

sc.pl.tsne(adata, color=['Sst', 'Gad1', 'Gad2', 'Npy', 'Tac1', 'Tac2', 'Htr3a', 'Slc6a4', 'Cck', 'Fos', 'Jun', 'Lhx6', 'Cacna2d2', 'Erbb4', 'Snhg11', 'Nxph2', 'Sox17', 'Arx', 'Pax6', 'Olig1', 'Meg3'], save="tSNE_adata_present_markers.png")

sc.pl.tsne(adata, color="louvain", save="adata_annotated.png", title="t-SNE CNS cell clusters from scRNA-seq", legend_loc="on data", legend_fontsize=6)

# sc.tl.umap(adata)
# with open("adata_var_names.out", 'w') as file:
#     for item in var_names1:
#         file.write(str(item) + "\n")

# sc.pp.calculate_qc_metrics(adata, expr_type='counts', var_type='genes', qc_vars=(), percent_top=(50, 100, 200, 500), inplace=True)

# # sc.pl.tsne(adata.uns, color='Pdyn', save="tSNE_pca2_Pdyn_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Sst', save="tSNE_pca2_Sst_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Oprk1', save="tSNE_pca2_Oprk1_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Slc17a6', save="tSNE_pca2_vGluT2_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Slc17a7', save="tSNE_pca2_vGluT1_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Slc17a8', save="tSNE_pca2_vGluT3_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Slc32a1', save="tSNE_pca2_vGluT2_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Htr3a', save="tSNE_pca2_Htr3a_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Drd1', save="tSNE_pca2_Drd1_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Drd2', save="tSNE_pca2_Drd2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Slc6a4', save="tSNE_pca2_Slc6a4_SERT_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Slc6a3', save="tSNE_pca2_Slc6a3_DAT_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Slc6a2', save="tSNE_pca2_Slc6a2_NET_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Crh', save="tSNE_pca2_Crh_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Cck', save="tSNE_pca2_Cck_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Chat', save="tSNE_pca2_ChAT_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Galp', save="tSNE_pca2_Galp_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Pomc', save="tSNE_pca2_Pomc_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Oxt', save="tSNE_pca2_Oxt_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Prkcd', save="tSNE_pca2_Prkcd_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Pnmt', save="tSNE_pca2_Pnmt_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Ddc', save="tSNE_pca2_Ddc_AADC_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Kcnmb1', save="tSNE_pca2_Kcnmb1_Slo1_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Kcnj3', save="tSNE_pca2_Kcnj3_GIRK1_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Kcnj6', save="tSNE_pca2_Kcnj6_GIRK2_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Kcnj9', save="tSNE_pca2_Kcnj9_GIRK3_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Kcnj5', save="tSNE_pca2_Kcnj5_GIRK4_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Tac1', save="tSNE_pca2_Tac1_SubP_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Tac2', save="tSNE_pca2_Tac2_neurokininB_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Avp', save="tSNE_pca2_Avp_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Tph1', save="tSNE_pca2_Tph1_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Tph2', save="tSNE_pca2_Tph2_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Avp', save="tSNE_pca2_Avp_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Fos', save="tSNE_pca2_Fos_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Jun', save="tSNE_pca2_Jun_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Npy', save="tSNE_pca2_Npy_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Vip', save="tSNE_pca2_Vip_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Esr1', save="tSNE_pca2_Esr1_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Penk', save="tSNE_pca2_Penk_louvain_adataMarkers.png")
# # sc.pl.tsne(adata, color='Dbh', save="tSNE_pca2_Dbh_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Gad2', save="tSNE_pca2_Gad2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Gad1', save="tSNE_pca2_Gad1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Gfap', save="tSNE_pca2_Gfap_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Gphn', save="tSNE_pca2_Gphn_Geph_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Dlg4', save="tSNE_pca2_Dlg4_Psd95_louvain_adata.png")
# sc.pl.tsne(adata, color='Rbfox3', save="tSNE_pca2_Rbfox3_NeuN_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Pvalb', save="tSNE_pca2_Pvalb_PV_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Syt2', save="tSNE_pca2_Syt2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Crhbp', save="tSNE_pca2_Crhbp_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Arx', save="tSNE_pca2_Arx_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Lhx6','Cacna2d2','Erbb4','4932438A13Rik','Snhg11','Plec','Cdc42ep2','Htt','App''Ranbp2','Tox3','Nrip3','Slc7a3','Npy2r','Nxph2','Nos1','Camk2n2','Adcy2','Kcnip1','Ank1','Sox6','Arx', save="tSNE_pca2_Cacna2d2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Sox6', save="tSNE_pca2_Sox6_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Ank1', save="tSNE_pca2_Ank1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Kcnip1', save="tSNE_pca2_Kcnip1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Adcy2', save="tSNE_pca2_Adcy2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Camk2n2', save="tSNE_pca2_Camk2m2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Nos1', save="tSNE_pca2_Nos1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Nxph2', save="tSNE_pca2_Nxph2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Npy2r', save="tSNE_pca2_Camk2m2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Slc7a3', save="tSNE_pca2_Slc7a3_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Nrip3', save="tSNE_pca2_Nrip3_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Tox3', save="tSNE_pca2_Tox3_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Ranbp2', save="tSNE_pca2_Ranbp2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Hivep2', save="tSNE_pca2_Hivep2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Htt', save="tSNE_pca2_Htt_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Cdc42ep2', save="tSNE_pca2_Cdc42ep2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Plec', save="tSNE_pca2_Plec_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Snhg11', save="tSNE_pca2_Snhg11_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='4932438A13Rik', save="tSNE_pca2_4932438A13Rik_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Erbb4', save="tSNE_pca2_Erbb4_louvain_adataMarkers.png")

# sc.tl.louvain(adata)
# sc.pl.tsne(adata, color='louvain', save="tSNE_pca2_louvain_adataMarkers.png")

# # sc.pl.umap(adata, color='louvain', save="uMAP_pca2_louvain_adata.png")
# # sc.pl.umap(adata, color='Sst', save="uMAP_pca2_louvain_adata.png")
# # sc.pl.umap(adata, color='Lhx6', save="uMAP_pca2_Lhx_louvain_adata.png")
# # sc.pl.umap(adata, color='Gad1', save="uMAP_pca2_Gad1_louvain_adata.png")
# # sc.pl.umap(adata, color='Gad2', save="uMAP_pca2_Gad2_louvain_adata.png")
# sc.tl.rank_genes_groups(adata, groupby='louvain', method="logreg")
# sc.pl.rank_genes_groups(adata, save="RankGenes_PCA50_LogReg_Sst_louvain_adata.png")
# sc.pl.rank_genes_groups_heatmap(adata, groupby=['louvian'], save="RankGenesPCA50_heatmap_LogReg_louvain_adata.png")
# # sc.tl.rank_genes_groups(adata, groupby=['louvain'], method="t-test")
# # sc.pl.rank_genes_groups(adata, save="RankGenes_ttest_louvain_adata1-3.png")
# # sc.pl.rank_genes_groups_heatmap(adata, groupby='louvian', save="RankGenes_heatmap_t-test_louvain_adata.png")
# # sc.pl.clustermap(adata, obs_key='louvain', save="clustermap_zheng17_louvain.png")
# # sc.pl.clustermap(adata, obs_key='louvain', save="clustermap_seurat_louvain.png")
# # print(var_names1)
# # HTR3A
# # CHRM3
# # NPY1R
# # VIPR1
# # RGS12
# # ADCY9
# # SYT10
# # PNOC
# # CRH
# # TAC2
# # HS6ST3
# # ADCY2
# # PDE2A "SST"

# # for i in filtering_method:
# # 	print(i) 
# # for i in var_names1:
# # 	print(i) 

# adata = sc.read_10x_h5(sys.argv[1])
# sc.tl.pca(adata, n_comps=50)
# sc.pp.recipe_seurat(adata)
# sc.pp.neighbors(adata)
# sc.tl.louvain(adata)


# # sc.tl.dotplot(adata)
# sc.tl.tsne(adata)
# # sc.tl.stacked_violin(adata)
# sc.tl.umap(adata)

# # sc.pl.tsne(adata, color='Pdyn', save="tSNE_pca2_Pdyn_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Sst', save="tSNE_pca2_Sst_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Oprk1', save="tSNE_pca2_Oprk1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Slc17a6', save="tSNE_pca2_vGluT2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Slc17a7', save="tSNE_pca2_vGluT1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Slc17a8', save="tSNE_pca2_vGluT3_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Slc32a1', save="tSNE_pca2_vGluT2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Htr3a', save="tSNE_pca2_Htr3a_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Drd1', save="tSNE_pca2_Drd1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Drd2', save="tSNE_pca2_Drd2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Slc6a4', save="tSNE_pca2_Slc6a4_SERT_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Slc6a3', save="tSNE_pca2_Slc6a3_DAT_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Slc6a2', save="tSNE_pca2_Slc6a2_NET_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Crh', save="tSNE_pca2_Crh_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Cck', save="tSNE_pca2_Cck_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Chat', save="tSNE_pca2_ChAT_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Galp', save="tSNE_pca2_Galp_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Pomc', save="tSNE_pca2_Pomc_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Oxt', save="tSNE_pca2_Oxt_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Prkcd', save="tSNE_pca2_Prkcd_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Pnmt', save="tSNE_pca2_Pnmt_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Ddc', save="tSNE_pca2_Ddc_AADC_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Kcnmb1', save="tSNE_pca2_Kcnmb1_Slo1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Kcnj3', save="tSNE_pca2_Kcnj3_GIRK1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Kcnj6', save="tSNE_pca2_Kcnj6_GIRK2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Kcnj9', save="tSNE_pca2_Kcnj9_GIRK3_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Kcnj5', save="tSNE_pca2_Kcnj5_GIRK4_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Tac1', save="tSNE_pca2_Tac1_SubP_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Tac2', save="tSNE_pca2_Tac2_neurokininB_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Avp', save="tSNE_pca2_Avp_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Tph1', save="tSNE_pca2_Tph1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Tph2', save="tSNE_pca2_Tph2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Avp', save="tSNE_pca2_Avp_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Fos', save="tSNE_pca2_Fos_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Jun', save="tSNE_pca2_Jun_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Npy', save="tSNE_pca2_Npy_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Vip', save="tSNE_pca2_Vip_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Esr1', save="tSNE_pca2_Esr1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Penk', save="tSNE_pca2_Penk_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Dbh', save="tSNE_pca2_Dbh_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Gad2', save="tSNE_pca2_Gad2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Gad1', save="tSNE_pca2_Gad1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Gad2', save="tSNE_pca2_Gad2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Gfap', save="tSNE_pca2_Gfap_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Gphn', save="tSNE_pca2_Gphn_Geph_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Dlg4', save="tSNE_pca2_Dlg4_Psd95_louvain_adata.png")
# sc.pl.tsne(adata, color='Rbfox3', save="tSNE_pca2_Rbfox3_NeuN_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Pvalb', save="tSNE_pca2_Pvalb_PV_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Syt2', save="tSNE_pca2_Syt2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Crhbp', save="tSNE_pca2_Crhbp_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Arx', save="tSNE_pca2_Arx_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Cacna2d2', save="tSNE_pca2_Cacna2d2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Sox6', save="tSNE_pca2_Sox6_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Ank1', save="tSNE_pca2_Ank1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Kcnip1', save="tSNE_pca2_Kcnip1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Adcy2', save="tSNE_pca2_Adcy2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Camk2n2', save="tSNE_pca2_Camk2m2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Nos1', save="tSNE_pca2_Nos1_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Nxph2', save="tSNE_pca2_Nxph2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Npy2r', save="tSNE_pca2_Camk2m2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Slc7a3', save="tSNE_pca2_Slc7a3_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Nrip3', save="tSNE_pca2_Nrip3_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Tox3', save="tSNE_pca2_Tox3_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Ranbp2', save="tSNE_pca2_Ranbp2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Hivep2', save="tSNE_pca2_Hivep2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Htt', save="tSNE_pca2_Htt_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Cdc42ep2', save="tSNE_pca2_Cdc42ep2_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Plec', save="tSNE_pca2_Plec_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Snhg11', save="tSNE_pca2_Snhg11_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='4932438A13Rik', save="tSNE_pca2_4932438A13Rik_louvain_adataMarkers.png")
# sc.pl.tsne(adata, color='Erbb4', save="tSNE_pca2_Erbb4_louvain_adataMarkers.png")

# sc.pl.umap(adata, color='louvain', save="uMAP_pca2zheng17_louvain_adata2.png")
# sc.pl.umap(adata, color='Sst', save="uMAP_pca2_zheng17_louvain_adata2.png")
# sc.pl.umap(adata, color='Lhx6', save="uMAP_pca2_Lhx6zheng17_louvain_adata2.png")
# sc.pl.umap(adata, color='Gad1', save="uMAP_pca2_Gad1_zheng17_louvain_adata2.png")
# sc.pl.umap(adata, color='Gad2', save="uMAP_pca2_Gad2_zheng17_louvain_adata2.png")
# sc.pl.umap(adata, color='Gad1', save="uMAP_pca2_Gad1_zheng17_louvain_adata2.png")
# sc.pl.umap(adata, color='Gad2', save="uMAP_pca2_Gad2_zheng17_louvain_adata2.png")
# sc.pl.umap(adata, color='Pdyn', save="uMAP_pca2_zheng17_louvain_adata2.png")

# # sc.pl.stacked_violin(adata, var_names=var_names1, color='louvain', save="stacked_violin_zheng17_louvain_ranked_adata.png" )
# # sc.pl.dotplot(adata, var_names=var_names1, color='louvain', save="dot_plot_zheng17_louvain_adata.png" )
# sc.pl.tsne(adata, color='louvain', save="tSNE_pca2_zheng17_louvain_adata-3.png")
# sc.pl.tsne(adata, color='Sst', save="tSNE_pca2_Sst_zheng17_louvain_adata-3.png")
# sc.pl.tsne(adata, color='Lhx6', save="tSNE_pca2_Lhx6_zheng17_louvain_adata-2.png")
# sc.pl.tsne(adata, color='Gad1', save="tSNE_pca2_Gad1_zheng17_louvain_adata-2.png")
# sc.pl.tsne(adata, color='Gad2', save="tSNE_pca2_Gad2_zheng17_louvain_adata.png")
# sc.pl.umap(adata, color='n_counts_all', save="uMAP_pca2_ncounts_all_zheng17_adata.png")
# sc.pl.umap(adata, color='n_counts', save="uMAP_pca2_ncounts_zheng17_adata.png")
# # sc.pl.tsne(adata, color='Pdyn', save="tSNE_Pdyn_zheng17_louvain_adata-2.png")
# # sc.pl.tsne(adata, color='Oprk1', save="tSNE_Oprk1_zheng17_louvain_adata-2.png")


# sc.tl.rank_genes_groups(adata, groupby='louvain', method="logreg")
# # sc.pl.rank_genes_groups(adata, save="RankGenes_LogReg_Sst_zheng17_louvain_adata2-3.png")
# # sc.pl.rank_genes_groups_heatmap(adata, groupby='louvian', save="RankGenes_heatmap_zheng17_LogReg_louvain_adata.png")
# # sc.pl.rank_genes_groups_stacked_violin(adata, var_names=var_names1, save="RankGenes_stackedviolin_zheng17_LogReg_louvain_adata2-3.png")
# # sc.pl.rank_genes_groups_dotplot(adata, var_names=var_names1, save="RankGenes_dotplot_zheng17_LogReg_louvain_adata2-3.png")
# # sc.pl.rank_genes_groups_matrixplot(adata, var_names=var_names1, save="RankGenes_matrix_zheng17_LogReg_louvain_adata2-3.png")

# sc.tl.rank_genes_groups(adata, groupby='louvain', method="t-test")
# # sc.pl.rank_genes_groups(adata, save="RankGenes_ttest_zheng17_louvain_adata1-3.png")
# # sc.pl.rank_genes_groups_heatmap(adata, groupby='louvian', save="RankGenes_heatmap_zheng17_t-test_louvain_adata.png")
# # sc.pl.rank_genes_groups_stacked_violin(adata, var_names=var_names1, save="RankGenes_zheng17_stackedviolin_t-test_louvain_adata2-3.png")
# # sc.pl.rank_genes_groups_dotplot(adata, var_names=var_names1, save="RankGenes_dotplot_zheng17_t-test_louvain_adata2-3.png")
# # sc.pl.rank_genes_groups_matrixplot(adata, var_names=var_names1, save="RankGenes_matrixzheng17__t-test_louvain_adata2-3.png")
# var_names1 = adata.var_names[0:]
# with open("adata_var_names_zheng17.out", 'w') as file:
#     for item in var_names1:
#         file.write(str(item) + "\n")
# # sc.tl.leiden(adata)
# # # sc.tl.dotplot(adata)
# # sc.tl.tsne(adata)
# # # sc.tl.stacked_violin(adata)

# # sc.tl.umap(adata)

# # sc.pl.umap(adata, color='leiden', save="uMAP_leiden_zheng17_adata1.png")
# # sc.pl.stacked_violin(adata, color='leiden', save="stacked_violin_zheng17_leiden_ranked_adata1.png" )
# # # sc.pl.dotplot(adata, color='leiden', save="dot_plot_zheng17_leiden_adata1.png" )
# # sc.pl.tsne(adata, color='leiden', save="tSNE_zheng17_leiden_adata1-2.png")
# # # sc.pl.tsne(adata, color='Slc17a6', save="tSNE_vGluT2_zheng17_leiden_adata-2.png")
# # # sc.pl.tsne(adata, color='Gad1', save="tSNE_Gad1_zheng17_leiden_adata-2.png")
# # # sc.pl.tsne(adata, color='Pdyn', save="tSNE_Pdyn_zheng17_leiden_adata-2.png")
# # # sc.pl.tsne(adata, color='Oprk1', save="tSNE_Oprk1_zheng17_leiden_adata-2.png")

# # sc.tl.rank_genes_groups(adata, groupby='leiden', method="logreg")
# # sc.pl.rank_genes_groups(adata, save="RankGenes_LogReg_leiden_zheng17adata1_2-3.png")
# # sc.pl.rank_genes_groups_heatmap(adata, var_names=var_names1, save="RankGenes_heatmap_zheng17_LogReg_leidendata1_2-3.png")
# # sc.pl.rank_genes_groups_stacked_violin(adata, var_names=var_names1, save="RankGenes_stackedviolin_zheng17_LogReg_leiden_adata1_2-3.png")
# # sc.pl.rank_genes_groups_dotplot(adata, var_names=var_names1, save="RankGenes_dotplot_zheng17_LogReg_leiden_adata1_2-3.png")
# # sc.pl.rank_genes_groups_matrixplot(adata, var_names=var_names1, save="RankGenes_matrixplot_zheng17_LogReg_leiden_adata1_2-3.png")

# # sc.tl.rank_genes_groups(adata, groupby='leiden', method="t-test")
# # sc.pl.rank_genes_groups(adata, save="RankGenes_ttest_leiden_zheng17_adata1_1-3.png")
# # sc.pl.rank_genes_groups_heatmap(adata, var_names=var_names1, save="RankGenes_heatmap_t-test_leiden_zheng17_adata1_2-3.png")
# # sc.pl.rank_genes_groups_stacked_violin(adata, var_names=var_names1, save="RankGenes_stackedviolin_zheng17_t-test_leiden_adata1_2-3.png")
# # sc.pl.rank_genes_groups_dotplot(adata, var_names=var_names1, save="RankGenes_dotplot_zheng17_t-test_leiden_adata1_2-3.png")
# # sc.pl.rank_genes_groups_matrixplot(adata, var_names=var_names1, save="RankGenes_matrix_zheng17_t-test_leiden_adata1_2-3.png")

# # adata1 = sc.read_10x_h5(sys.argv[1])
# # adata1.var_names_make_unique()

# # sc.tl.pca(adata1, n_comps=50)
# # sc.pp.recipe_weinreb17(adata1)
# # sc.pp.neighbors(adata1)

# # sc.tl.louvain(adata1)
# # # sc.tl.dotplot(adata1)
# # sc.tl.tsne(adata1)
# # # sc.tl.stacked_violin(adata1)

# # sc.tl.umap(adata1)

# # sc.pl.umap(adata1, color='louvain', save="uMAP_louvain_weinreb17_adata1.png")
# # sc.pl.stacked_violin(adata1, var_names=var_names1, color='louvain', save="stacked_violin_weinreb17_louvain_ranked_adata1.png" )
# # # sc.pl.dotplot(adata1, color='louvain', save="dot_plot_weinreb17_louvain_adata1.png" )
# # sc.pl.tsne(adata1, color='louvain', save="tSNE_weinreb17_louvain_adata1-2.png")
# # # sc.pl.tsne(adata1, color='Slc17a6', save="tSNE_vGluT2_zheng17_louvain_adata1-2.png")
# # # sc.pl.tsne(adata1, color='Gad1', save="tSNE_Gad1_weinreb17_louvain_adata1-2.png")
# # # sc.pl.tsne(adata1, color='Pdyn', save="tSNE_Pdyn_weinreb17_louvain_adata-2.png")
# # # sc.pl.tsne(adata1, color='Oprk1', save="tSNE_Oprk1_weinreb17_louvain_adata1-2.png")

# # sc.tl.rank_genes_groups(adata1, groupby='louvain', method="logreg")
# # sc.pl.rank_genes_groups(adata1, save="RankGenes_LogReg_louvain__weinreb17adata1_2-3.png")
# # # sc.pl.rank_genes_groups_heatmap(adata1, var_names=var_names1, save="RankGenes_heatmap_weinreb17_LogReg_louvain_adata1_2-3.png")
# # # sc.pl.rank_genes_groups_stacked_violin(adata1, var_names=var_names1, save="RankGenes_stackedviolin_weinreb17_LogReg_louvain_adata1_2-3.png")
# # # sc.pl.rank_genes_groups_dotplot(adata1, var_names=var_names1, save="RankGenes_dotplot_weinreb17_LogReg_louvain_adata1_2-3.png")
# # # sc.pl.rank_genes_groups_matrixplot(adata1, var_names=var_names1, save="RankGenes_matrixplot_weinreb17_LogReg_louvain_adata1_2-3.png")

# # sc.tl.rank_genes_groups(adata1, groupby='louvain', method="t-test")
# # sc.pl.rank_genes_groups(adata1, save="RankGenes_ttest_louvain_weinreb17_adata1_1-3.png")
# # sc.pl.rank_genes_groups_heatmap(adata1, var_names=var_names1, save="RankGenes_heatmap_t-test_louvain_weinreb17_adata1_2-3.png")
# # sc.pl.rank_genes_groups_stacked_violin(adata1, var_names=var_names1, save="RankGenes_stackedviolin_weinreb17_t-test_louvain_adata1_2-3.png")
# # sc.pl.rank_genes_groups_dotplot(adata1, var_names=var_names1, save="RankGenes_dotplot_weinreb17_t-test_louvain_adata1_2-3.png")
# # sc.pl.rank_genes_groups_matrixplot(adata1, var_names=var_names1, save="RankGenes_matrix_weinreb17_t-test_louvain_adata1_2-3.png")


# # sc.tl.leiden(adata1)
# # # sc.tl.dotplot(adata1)
# # sc.tl.tsne(adata1)
# # # sc.tl.stacked_violin(adata1)

# # sc.tl.umap(adata1)

# # sc.pl.umap(adata1, color='leiden', save="uMAP_leiden_weinreb17_adata1.png")
# # sc.pl.stacked_violin(adata1, var_names=var_names1, color='leiden', save="stacked_violin_weinreb17_leiden_ranked_adata1.png" )
# # # sc.pl.dotplot(adata1, var_names=var_names1, color='leiden', save="dot_plot_weinreb17_leiden_adata1.png" )
# # sc.pl.tsne(adata1, color='leiden', save="tSNE_weinreb17_leiden_adata1-2.png")
# # # sc.pl.tsne(adata1, color='Slc17a6', save="tSNE_vGluT2_weinreb17_leiden_adata1-2.png")
# # # sc.pl.tsne(adata1, color='Gad1', save="tSNE_Gad1_weinreb17_leiden_adata1-2.png")
# # # sc.pl.tsne(adata1, color='Pdyn', save="tSNE_Pdyn_weinreb17_leiden_adata1-2.png")
# # # sc.pl.tsne(adata1, color='Oprk1', save="tSNE_Oprk1_weinreb17_leiden_adata1-2.png")

# # sc.tl.rank_genes_groups(adata1, groupby='leiden', method="logreg")
# # sc.pl.rank_genes_groups(adata1, save="RankGenes_LogReg_leiden__weinreb17adata1_2-3.png")
# # # sc.pl.rank_genes_groups_heatmap(adata1, var_names=var_names1, save="RankGenes_heatmap_weinreb17_LogReg_leidendata1_2-3.png")
# # # sc.pl.rank_genes_groups_stacked_violin(adata1, var_names=var_names1, save="RankGenes_stackedviolin_weinreb17_LogReg_leiden_adata1_2-3.png")
# # # sc.pl.rank_genes_groups_dotplot(adata1, var_names=var_names1, save="RankGenes_dotplot_weinreb17_LogReg_leiden_adata1_2-3.png")
# # # sc.pl.rank_genes_groups_matrixplot(adata1, var_names=var_names1, save="RankGenes_matrixplot_weinreb17_LogReg_leiden_adata1_2-3.png")

# # sc.tl.rank_genes_groups(adata1, groupby='leiden', method="t-test")
# # sc.pl.rank_genes_groups(adata1, save="RankGenes_ttest_louvain_weinreb17_adata1_1-3.png")
# # sc.pl.rank_genes_groups_heatmap(adata1, var_names=var_names1,  save="RankGenes_heatmap_t-test_leiden_weinreb17_adata1_2-3.png")
# # sc.pl.rank_genes_groups_stacked_violin(adata1, var_names=var_names1, save="RankGenes_stackedviolin_weinreb17_t-test_leiden_adata1_2-3.png")
# # sc.pl.rank_genes_groups_dotplot(adata1, var_names=var_names1, save="RankGenes_dotplot_weinreb17_t-test_leiden_adata1_2-3.png")
# # sc.pl.rank_genes_groups_matrixplot(adata1, var_names=var_names1, save="RankGenes_matrix_weinreb17_t-test_leiden_adata1_2-3.png")

# adata2 = sc.read_10x_h5(sys.argv[1])
# adata2.var_names_make_unique()

# sc.tl.pca(adata2, n_comps=2)
# sc.pp.recipe_seurat(adata2)
# sc.pp.neighbors(adata2)

# sc.tl.louvain(adata2)
# # sc.tl.dotplot(adata2)
# sc.tl.tsne(adata2)
# # sc.tl.stacked_violin(adata2)

# sc.tl.umap(adata2)

# sc.pl.umap(adata2, color='n_counts_all', save="uMAP_pca2_ncounts_all_seurat_adata3.png")
# sc.pl.umap(adata2, color='n_counts', save="uMAP_pca2_ncounts_seurat_adata3.png")
# sc.pl.umap(adata2, color='louvain', save="uMAP_pca2_louvain_seurat_adata3.png")
# sc.pl.umap(adata2, color='Sst', save="uMAP_pca2_Sst_louvain_seurat_adata3.png")
# # sc.pl.umap(adata2, color='Gad1', save="uMAP_gad1_louvain_seurat_adata3.png")
# # sc.pl.umap(adata2, color='Gad2', save="uMAP_gad2_louvain_seurat_adata3.png")
# # sc.pl.umap(adata2, color='Pdyn', save="uMAP_Pdyn_louvain_seurat_adata3.png")
# sc.pl.umap(adata2, color='Lhx6', save="uMAP_pca2_Lhx6_louvain_seurat_adata2-4.png")
# # sc.pl.stacked_violin(adata2, var_names=var_names1,color='louvain', save="stacked_violin_weinreb17seurat_louvain_ranked_adata2.png" )
# # sc.pl.dotplot(adata2, var_names=var_names1, color='louvain', save="dot_plot_weinreb17seurat_louvain_adata2.png" )
# sc.pl.tsne(adata2, color='louvain', save="tSNE_pca2_seurat_louvain_adata2-4.png")
# sc.pl.tsne(adata2, color='Sst', save="tSNE_pca2_Sst_seurat_louvain_adata2-4.png")
# # sc.pl.tsne(adata2, color='Slc17a6', save="tSNE_vGluT2_seurat_louvain_adata2-2.png")
# # sc.pl.tsne(adata2, color='Gad1', save="tSNE_Gad1_seurat_louvain_adata2-3.png")
# # sc.pl.tsne(adata2, color='Gad2', save="tSNE_Gad1_seurat_louvain_adata2-3.png")
# # sc.pl.tsne(adata2, color='Pdyn', save="tSNE_Pdyn_seurat_louvain_adata2-3.png")
# sc.pl.tsne(adata2, color='Lhx6', save="tSNE_pca2_Lhx6_seurat_louvain_adata2-3.png")
# # sc.pl.tsne(adata2, color='Oprk1', save="tSNE_Oprk1_seurat_louvain_adata2-2.png")


# # sc.tl.rank_genes_groups(adata2, groupby='louvain', method="logreg")
# # sc.pl.rank_genes_groups(adata2, save="RankGenes_LogReg_louvain_seuratadata2_2-3.png")
# # sc.pl.umap(adata2, color='Gad1', save="uMAP_gad1_louvain_seurat_adata2-4.png")
# # sc.pl.umap(adata2, color='Gad2', save="uMAP_gad2_louvain_seurat_adata2-4.png")
# # sc.pl.rank_genes_groups_heatmap(adata2, groupby='louvian', save="RankGenes_heatmap_weinreb17seurat_LogReg_louvain_adata2.png")
# # sc.pl.rank_genes_groups_stacked_violin(adata2, var_names=var_names1, save="RankGenes_stackedviolin_weinreb17seurat_LogReg_louvain_adata2_2-3.png")
# # sc.pl.rank_genes_groups_dotplot(adata2, var_names=var_names1,save="RankGenes_dotplot_weinreb17seurat_LogReg_louvain_adata2_2-3.png")
# # sc.pl.rank_genes_groups_matrixplot(adata2, var_names=var_names1, save="RankGenes_matrixplot_seurat_LogReg_louvain_adata2_2-3.png")

# # sc.tl.rank_genes_groups(adata2, groupby='louvain', method="t-test")
# # sc.pl.rank_genes_groups(adata2, save="RankGenes_ttest_louvain_seurat_adata2_2-3.png")
# # sc.pl.rank_genes_groups_heatmap(adata2, groupby='louvian', save="RankGenes_heatmap_t-test_louvain_seurat_adata2.png")
# # sc.pl.rank_genes_groups_stacked_violin(adata2, var_names=var_names1, save="RankGenes_stackedviolin_seurat_t-test_louvain_adata2_2-3.png")
# # sc.pl.rank_genes_groups_dotplot(adata2, var_names=var_names1, save="RankGenes_dotplot_seurat_t-test_louvain_adata2_2-3.png")
# # sc.pl.rank_genes_groups_matrixplot(adata2, var_names=var_names1, save="RankGenes_matrix_seurat_t-test_louvain_adata2_2-3.png")
# sc.pl.umap(adata2, color='Gad1', save="uMAP_gad1_louvain_seurat_adata2-4.png")
# sc.pl.umap(adata2, color='Gad2', save="uMAP_gad2_louvain_seurat_adata2-4.png")

# with open("adata_var_names_seurat.out", 'w') as file:
#     for item in var_names1:
#         file.write(str(item) + "\n")
# # sc.tl.leiden(adata2)
# # # sc.tl.dotplot(adata2)
# # sc.tl.tsne(adata2)
# # # sc.tl.stacked_violin(adata2)

# # sc.tl.umap(adata2)

# # sc.pl.umap(adata2, color='leiden', save="uMAP_leiden_seurat_adata2.png")
# # # sc.pl.stacked_violin(adata2, color='leiden', save="stacked_violin_seurat_leiden_ranked_adata2.png" )
# # # sc.pl.dotplot(adata2, color='leiden', save="dot_plot_seurat_leiden_adata2.png" )
# # sc.pl.tsne(adata2, color='leiden', save="tSNE_seurat_leiden_adata2-2.png")
# # # sc.pl.tsne(adata2, color='Slc17a6', save="tSNE_vGluT2_seurat_leiden_adata2-2.png")
# # # sc.pl.tsne(adata2, color='Gad1', save="tSNE_Gad1_seurat_leiden_adata2-2.png")
# # # sc.pl.tsne(adata2, color='Pdyn', save="tSNE_Pdyn_seurat_leiden_adata2-2.png")
# # # sc.pl.tsne(adata2, color='Oprk1', save="tSNE_Oprk1_seurat_leiden_adata2-2.png")

# # sc.tl.rank_genes_groups(adata2, groupby='leiden', method="logreg")
# # sc.pl.rank_genes_groups(adata2, save="RankGenes_LogReg_leiden__seuratadata2_2-3.png")
# # # sc.pl.rank_genes_groups_heatmap(adata2, var_names=var_names1, save="RankGenes_heatmap_seurat_LogReg_leidendata2_2-3.png")
# # # sc.pl.rank_genes_groups_stacked_violin(adata2, var_names=var_names1, save="RankGenes_stackedviolin_seurat_LogReg_leiden_adata2_2-3.png")
# # # sc.pl.rank_genes_groups_dotplot(adata2, var_names=var_names1, save="RankGenes_dotplot_seurat_LogReg_leiden_adata2_2-3.png")
# # # sc.pl.rank_genes_groups_matrixplot(adata2, var_names=var_names1, save="RankGenes_matrixplot_seurat_LogReg_leiden_adata2_2-3.png")

# # sc.tl.rank_genes_groups(adata2, groupby='leiden', method="t-test")
# # sc.pl.rank_genes_groups(adata2, save="RankGenes_ttest_leiden_seurat_adata2_1-3.png")
# # # sc.pl.rank_genes_groups_heatmap(adata2, var_names=var_names1, save="RankGenes_heatmap_t-test_leiden_seurat_adata2_2-3.png")
# # sc.pl.rank_genes_groups_stacked_violin(adata2, var_names=var_names1 ,groupby='leiden' , save="RankGenes_stackedviolin_seurat_t-test_leiden_adata2_2-3.png")
# # sc.pl.rank_genes_groups_dotplot(adata2, var_names=var_names1, save="RankGenes_dotplot_seurat_t-test_leiden_adata2_2-3.png")
# # sc.pl.rank_genes_groups_matrixplot(adata2, var_names=var_names1, save="RankGenes_matrix_seurat_t-test_leiden_adata2_2-3.png")
# # sc.pl.rank_genes_groups(adata, save="RankGenes_LogReg_louvain_adata2.png")

# # sc.tl.umap(adata)
# # sc.pl.umap(adata, color='louvain' , save="uMAP_zheng17_louvain_adata.png")

# # sc.tl.tsne(adata)
# # sc.pl.tsne(adata, color='rank_genes_groups', save="tSNE_zheng17_louvain_adata-1.png")


# # adata1 = sc.tl.rank_genes_groups(adata, groupby='louvain', method="t-test")
# # sc.pl.rank_genes_groups(adata1, save="RankGenes_ttest_louvain_adata1.png")

# # adata2 = sc.tl.rank_genes_groups(adata, groupby='louvain', method="logreg")
# # sc.pl.rank_genes_groups(adata1, save="RankGenes_LogReg_louvain_adata2.png")

# # sc.tl.umap(adata1)
# # sc.pl.umap(adata1, color='louvain', save="uMAP_zheng17_louvain_adata.png")

# # sc.tl.tsne(adata1)
# # sc.pl.tsne(adata1, color='louvain', save="tSNE_zheng17_louvain_adata.png")

# # sc.tl.umap(adata2)
# # sc.pl.umap(adata2, save="uMAP_zheng17_louvain_adata.png")

# # sc.tl.tsne(adata2)
# # sc.pl.tsne(adata2, color='louvain', save="tSNE_zheng17_louvain_adata.png")
# # sc.pl.dotplot(adata, var_names=var_names1, color='louvain', save="dot_plot_zheng17_louvain_adata.png" )
# # sc.tl.rank_genes_groups(adata, groupby='Sst', method="logreg")

# sc.pl.clustermap(adata, obs_key='louvain', save="clustermap_zheng17_louvain.png")
# sc.pl.clustermap(adata, obs_key='louvain', save="clustermap_seurat_louvain.png")



