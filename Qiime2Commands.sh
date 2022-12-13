#We can take our output from SBAnalyzer and move into Qiime2

#Activate qiime2
conda activate qiime2-2021.8

#Next, import the data into qiime2
qiime tools import \
  --input-path sequence_table_merged.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path table.qza

qiime tools import \
  --input-path ASVs_merged.fasta \
  --output-path rep-seqs.qza \
  --type 'FeatureData[Sequence]'
  
qiime tools import \
  --input-path taxonomy_merged.tax \
  --output-path taxonomy.qza \
  --type 'FeatureData[Taxonomy]'

#Create phylogenetic trees based on the rep-seqs file.
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

#Filter the feature table to include only the samples of interest by using a metadata file that only contains the samples we want.
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metadata.txt \
  --o-filtered-table filtered-table.qza

#Create a visualization of the table
qiime feature-table summarize \
  --i-table filtered-table.qza \
  --o-visualization filtered-table.qzv

#Calculate alpha and beta diversity metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table filtered-table.qza \
  --p-sampling-depth 3746 \
  --m-metadata-file metadata.txt \
  --output-dir core-metrics-results

#Calculate statistics for alpha diversity based on metadata columns
#Faith's phylogenetic diversity
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file metadata.txt \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv	

#Pielou's evenness
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file metadata.txt \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

#By default, the core-metrics function doesn't include Shannon diversity. We'll calculate that separately.
qiime diversity alpha \
  --i-table core-metrics-results/rarefied_table.qza \
  --p-metric shannon \
  --o-alpha-diversity core-metrics-results/shannon.qza

#Now, calculate the statistics for Shannon diversity based on metadata columns.
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon.qza \
  --m-metadata-file metadata.txt \
  --o-visualization core-metrics-results/shannon.qzv


#Calculate statistics for phylogeny-based beta diversity based on the "Country" metadata column

#Unweighted UniFrac Permanova by country
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column Country \
  --o-visualization core-metrics-results/unweighted-unifrac-Country.qzv \
  --p-pairwise

#Unweighted UniFrac dispersion test (PERMDISP) by country
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column Country \
  --o-visualization core-metrics-results/permdisp-unweighted-unifrac-Country.qzv \
  --p-method permdisp \
  --p-pairwise

#Weighted UniFrac Permanova by country
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column Country \
  --o-visualization core-metrics-results/weighted-unifrac-Country.qzv \
  --p-pairwise

#Weighted UniFrac dispersion test (PERMDISP) by country
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column Country \
  --o-visualization core-metrics-results/permdisp-weighted-unifrac-Country.qzv \
  --p-method permdisp \
  --p-pairwise

#Calculate statistics for the non-phylogenetic beta diversity based on the "Country" metadata column 

#Due to a batching effect from DADA2, Jaccard and Bray-Curtis show artifacts. To solve this, collapse to the strain level.
mkdir TaxaLevel8
qiime taxa collapse \
 --i-table filtered-table.qza \
 --i-taxonomy taxonomy.qza \
 --p-level 8 \
 --o-collapsed-table table-8.qza

#Create the core metrics directory with the alpha and beta diversity calculations using the collapsed table.
qiime diversity core-metrics \
  --i-table table-8.qza \
  --p-sampling-depth 3746 \
  --m-metadata-file metadata.txt \
  --output-dir TaxaLevel8/core-metrics-results

#Bray-Curtis Permanova by country
qiime diversity beta-group-significance \
  --i-distance-matrix TaxaLevel8/core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column Country \
  --o-visualization TaxaLevel8/core-metrics-results/bray-curtis-Country.qzv \
  --p-pairwise

#Jaccard Permanova by country
qiime diversity beta-group-significance \
  --i-distance-matrix TaxaLevel8/core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column Country \
  --o-visualization TaxaLevel8/core-metrics-results/jaccard-Country.qzv \
  --p-pairwise

#Bray-Curtis dispersion test (PERMDISP) by country
qiime diversity beta-group-significance \
  --i-distance-matrix TaxaLevel8/core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column Country \
  --o-visualization TaxaLevel8/core-metrics-results/bray-curtis-Country-permdisp.qzv \
  --p-method permdisp \
  --p-pairwise

#Jaccard dispersion test (PERMDISP) by country
qiime diversity beta-group-significance \
  --i-distance-matrix TaxaLevel8/core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column Country \
  --o-visualization TaxaLevel8/core-metrics-results/jaccard-Country-permdisp.qzv \
  --p-method permdisp \
  --p-pairwise

#Next, we will use ANCOM to calculate taxa that differ significantly between the US and Colombia.
#We will do this from the phylum level to the strain level
taxlevels="2 3 4 5 6 7 8"

#Make a directory for the ANCOM comparisons, and copy the necessary files.
mkdir ANCOM
mkdir ANCOM/Country
cp filtered-table.qza ANCOM/Country
cp metadata.txt ANCOM/Country
cp taxonomy.qza ANCOM/Country
cd ANCOM/Country

		
for level in $taxlevels
do
	qiime taxa collapse \
		--i-table filtered-table.qza \
 		--i-taxonomy taxonomy.qza \
 		--p-level ${level} \
  		--o-collapsed-table table-${level}.qza

	qiime composition add-pseudocount \
		--i-table table-${level}.qza \
		--o-composition-table comp-table-${level}.qza
  
  	qiime composition ancom \
 		--i-table comp-table-${level}.qza \
		--m-metadata-file metadata.txt \
		--m-metadata-column Country \
		--o-visualization ${level}-ancom.qzv
done
cd ../..

#We can now move to R to generate visualizations.