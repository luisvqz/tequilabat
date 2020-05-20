#Tmux is recomended for long run in a remote server
tmux
#neccesary files:
#data: a single folder with all the fastq files
#mapping file: a single folder with one file mapfile20200504.txt
#taxo: a single folder with the a trained classifier silva-132-99-515-806-nb-classifier.qza
cd Documents/Migration/
ls
source activate qiime2-2020.2

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path data \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path outputff/demux-paired-end.qza 

qiime demux summarize \
--i-data outputff/demux-paired-end.qza \
--o-visualization outputff/demux-paired-end.qzv

qiime tools view outputff/demux-paired-end.qzv

qiime dada2 denoise-paired \
--i-demultiplexed-seqs outputff/demux-paired-end.qza \
--p-trim-left-f 23 \
--p-trim-left-r 20 \
--p-trunc-len-f 200 \
--p-trunc-len-r 200 \
--p-n-threads 20 \
--output-dir outputff/DADA2stats \
--o-representative-sequences outputff/rep-seqs-dada2.qza \
--o-table outputff/table-dada2.qza \
--verbose

# denoising results 

qiime feature-table summarize \
--i-table outputff/table-dada2.qza \
--o-visualization outputff/table-dada2.qzv \
--m-sample-metadata-file mapping/Mapfile20200504.txt

#denoising summary by step

qiime metadata tabulate \
  --m-input-file outputff/DADA2stats/denoising_stats.qza \
  --o-visualization outputff/DADA2stats/denoising_stats.qzv

#ASVs in dataset

qiime feature-table tabulate-seqs \
  --i-data outputff/rep-seqs-dada2.qza \
  --o-visualization outputff/rep-seqs-dada2.qzv

#Taxonomy assigment using Silva 132-99-515-816

#####SILVA
qiime feature-classifier classify-sklearn \
  --i-classifier taxo/silva-132-99-515-806-nb-classifier.qza \
  --i-reads outputff/rep-seqs-dada2.qza \
  --p-reads-per-batch 1000 \
  --o-classification outputff/taxonomy.qza \
  --verbose

## Taxomony as qzv

qiime metadata tabulate \
  --m-input-file outputff/taxonomy.qza \
  --o-visualization outputff/taxonomy.qzv

#Check results in boxplot (prior to plastid and mitochondria filtering)

qiime taxa barplot \
  --i-table outputff/table-dada2.qza \
  --i-taxonomy outputff/taxonomy.qza \
  --m-metadata-file mapping/Mapfile20200504.txt \
  --o-visualization outputff/taxa-bar-plotssectionunf.qzv

### filtering mitochodria and chloroplasts
#for dada2 table
qiime taxa filter-table \
  --i-table outputff/table-dada2.qza  \
  --i-taxonomy outputff/taxonomy.qza \
  --p-exclude Archaea,Mitochondra,Chloroplast,Unassigned \
  --p-include D_ \
  --o-filtered-table outputff/table-filtered0.qza

qiime taxa filter-table \
  --i-table outputff/table-filtered0.qza  \
  --i-taxonomy outputff/taxonomy.qza \
  --p-exclude D_4__Mitochondria \
  --p-include D_ \
  --o-filtered-table outputff/table-filtered.qza

#for representative sequences
qiime taxa filter-seqs \
  --i-sequences outputff/rep-seqs-dada2.qza \
  --i-taxonomy outputff/taxonomy.qza \
  --p-exclude Archaea,Mitochondra,Chloroplast,Unassigned \
  --p-include D_ \
  --o-filtered-sequences outputff/rep-seqs-filtered0.qza

qiime taxa filter-seqs \
  --i-sequences outputff/rep-seqs-filtered0.qza \
  --i-taxonomy outputff/taxonomy.qza \
  --p-exclude D_4__Mitochondria \
  --p-include D_ \
  --o-filtered-sequences outputff/rep-seqs-filtered.qza

#boxplots after filtering

qiime feature-table summarize \
--i-table outputff/table-filtered.qza \
--o-visualization outputff/table-filtered.qzv \
--m-sample-metadata-file mapping/Mapfile20200504.txt

#representative sequences after filtering

qiime feature-table tabulate-seqs \
  --i-data outputff/rep-seqs-filtered.qza \
  --o-visualization outputff/rep-seqs-filtered.qzv

#######tree for filtered ASVs
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences outputff/rep-seqs-filtered.qza \
  --o-alignment outputff/aligned-rep-seqs-filtered.qza \
  --o-masked-alignment outputff/masked-aligned-rep-seqs-filtered.qza \
  --o-tree outputff/unrooted-treefil.qza \
  --o-rooted-tree outputff/rooted-treefil.qza

qiime tools export \
  --input-path rooted-treefil.qza \
  --output-path exported-tree

#bloxplots for filtered

qiime taxa barplot \
  --i-table outputff/table-filtered.qza \
  --i-taxonomy outputff/taxonomy.qza \
  --m-metadata-file mapping/Mapfile20200504.txt \
  --o-visualization outputff/filtered-taxa-bar-plots.qzv

qiime tools view outputff/filtered-taxa-bar-plots.qzv

#Ok now we will export the data to be analysed into another sofware into R:

qiime tools export \
  --input-path outputff/table-filtered.qza \
  --output-path outputff/exported-table

qiime tools export \
  --input-path outputff/taxonomy.qza \
  --output-path outputff/taxonomy

biom convert -i outputff/exported-table/feature-table.biom -o outputff/exported-table/migration.txt --header-key taxonomy --to-tsv

#movimg into R
R
getwd()
table<-read.csv("exported-table/migration.txt",sep='\t',check.names=FALSE,skip=1)
head(table)
names(table)
Taxonomy<-read.csv("taxonomy/taxonomy.csv", header=TRUE)
head(Taxonomy)
names(Taxonomy)
table$taxonomy<-with(Taxonomy,Taxon[match(table$"#OTU ID",Taxonomy$Feature.ID)])
write.table(table,"exported-table/leptoTaxonomy.txt",row.names=FALSE,sep="\t")
quit()
n

cd exported-table
sed -i '1s/^/# Constructed from biom file\n/' leptoTaxonomy.txt
sed -i -e 's/"//g' leptoTaxonomy.txt
ls
subl exported-table/leptoTaxonomy.txt
cd ..
biom convert -i exported-table/leptoTaxonomy.txt -o exported-table/leptoTaxonomy.biom --table-type "OTU table" --process-obs-metadata taxonomy --to-json