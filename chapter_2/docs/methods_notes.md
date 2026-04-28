# Get genome data:
All genome data is access from: http://aequ.reefgenomics.org/download/

wget http://aequ.reefgenomics.org/download/equina_smart.rnam-trna.merged.ggf.curated.remredun.proteins.gff3.gz

## gunzip
gunzip equina_smart.rnam-trna.merged.ggf.curated.remredun.proteins.gff3.gz

## move to new file, better name
cp equina_smart.rnam-trna.merged.ggf.curated.remredun.proteins.gff3  a_equina.gff3
