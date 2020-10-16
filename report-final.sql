############################################################
# The sql script to organize tables for the report 
# and link gene level external /projects/analysis/report/annotations
###########################################################
# Variant impacts are decribed here
# https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
############################################################

## Concatenate all variant effects separately for refseq and for ensembl
DROP TABLE IF EXISTS grouped_impacts_refseq;
CREATE TABLE grouped_impacts_refseq AS select variant_id, COUNT(variant_id) AS Variant_impacts_num_refseq, 
  GROUP_CONCAT( distinct ( '[' || gene || ';' || impact_so || ';' || impact_severity || ';' || hgvsc || ';' || hgvsp || ']')) AS Grouped_variant_impacts_refseq
    FROM variant_impacts WHERE NOT INSTR(transcript,"EN") GROUP BY variant_id ;

DROP TABLE IF EXISTS grouped_impacts_enst;
CREATE TABLE grouped_impacts_enst AS select variant_id, COUNT(variant_id) AS Variant_impacts_num_enst, 
  GROUP_CONCAT( distinct ( '[' || gene || ';' || impact_so || ';' || impact_severity || ';' || hgvsc || ';' || hgvsp || ']')) AS Grouped_variant_impacts_enst 
    FROM variant_impacts WHERE INSTR(transcript,"EN") GROUP BY variant_id;

# distinct variant identifiers
DROP TABLE IF EXISTS var_ids;
CREATE TABLE var_ids AS SELECT DISTINCT variant_id
    FROM variant_impacts;

# Join grouped impacts
DROP TABLE IF EXISTS grouped_impacts_r;
CREATE TABLE grouped_impacts_r AS SELECT v.variant_id as variant_id,
  r.Variant_impacts_num_refseq,
  r.Grouped_variant_impacts_refseq
FROM var_ids v LEFT JOIN grouped_impacts_refseq r ON v.variant_id=r.variant_id;

DROP TABLE IF EXISTS grouped_impacts;
CREATE TABLE grouped_impacts AS SELECT v.variant_id as variant_id,
  v.Variant_impacts_num_refseq,
  v.Grouped_variant_impacts_refseq,
  r.Variant_impacts_num_enst,
  r.Grouped_variant_impacts_enst
FROM grouped_impacts_r v LEFT JOIN grouped_impacts_enst r ON v.variant_id=r.variant_id;

# Join genotypes
DROP TABLE IF EXISTS grouped_impacts_gts;
CREATE TABLE grouped_impacts_gts AS SELECT g.*,
  v.Variant_impacts_num_refseq,
  v.Grouped_variant_impacts_refseq,
  v.Variant_impacts_num_enst,
  v.Grouped_variant_impacts_enst
FROM genotypes g, grouped_impacts v where g.variant_id=v.variant_id;


## Select fields into the report
## Select fields from grouped variant effects : refseq and ensembl

DROP TABLE IF EXISTS variant_fields;
CREATE TABLE variant_fields AS SELECT 
        v.hgvsg AS Hgvsg,       
        CASE SUBSTR(r.gts,1,INSTR(r.gts,'/')-1)=SUBSTR(r.gts,INSTR(r.gts,'/')+1)
          WHEN 1 THEN 'hom'
          WHEN 0 THEN 'het'
        END AS Zygocity,

        r.gts AS Genotype,
        r.gt_alt_freqs AS Alt_Variant_Freq,
        v.dp AS Read_Depth,
        r.gt_alt_depths AS Alt_Read_Depth,
        r.gt_ref_depths || ',' || r.gt_alt_depths AS Allelic_Depths,

        v.variant_id AS variant_id,
        v.gene AS Gene,
        v.ref || '>'|| r.gts AS Variant,

        v.chrom AS Chrom,
        v.start+1 AS Coordinate,
        v.ref AS Ref,
        v.alt AS Alt,
        v.type AS Type,
        v.sub_type AS Sub_Type,
        
        CASE v.is_exonic
          WHEN 1 THEN 'yes'
          WHEN 0 THEN 'no'
        END AS Exonic,
        REPLACE(v.exon,'/','_') AS Exon,
        REPLACE(v.intron,'/','_') AS Intron,

        v.num_hom_ref AS Num_hom_ref,
        v.num_het AS Num_het,
        v.num_hom_alt AS Num_hom_alt,
        v.num_unknown AS Num_unknown,

        CAST(v.qual AS REAL) AS Quality,
        v.hotc_hgmd AS HGMD,
        v.hotc_localvariant AS Localvariant,
        SUBSTR(v.hotc_localvariant,INSTR(v.hotc_localvariant,"Count:")+LENGTH("Count:"),INSTR(v.hotc_localvariant,"|Freq:")-(INSTR(v.hotc_localvariant,"Count:")+LENGTH("Count:"))) AS Localvariant_Count,
        SUBSTR(v.hotc_localvariant,INSTR(v.hotc_localvariant,"Freq:")+LENGTH("Freq:"),INSTR(v.hotc_localvariant,"|Zyg:")-(INSTR(v.hotc_localvariant,"Freq:")+LENGTH("Freq:"))) AS Localvariant_Freq,
        SUBSTR(v.hotc_gl,INSTR(v.hotc_gl,"ENTREZ:")+LENGTH("ENTREZ:"),INSTR(v.hotc_gl,"|NM:")-(INSTR(v.hotc_gl,"ENTREZ:")+LENGTH("ENTREZ:"))) AS ENTREZ,
        '=HYPERLINK("https://www.ncbi.nlm.nih.gov/gene/?term=' || SUBSTR(v.hotc_gl,INSTR(v.hotc_gl,"ENTREZ:")+LENGTH("ENTREZ:"),INSTR(v.hotc_gl,"|NM:")-(INSTR(v.hotc_gl,"ENTREZ:")+LENGTH("ENTREZ:"))) || '")' AS linkNCBI,
        v.ensembl_gene_id AS VEPGeneID,
        SUBSTR(v.hotc_gl,INSTR(v.hotc_gl,"Description:")+LENGTH("Description:"),INSTR(v.hotc_gl,"|MIM:")-(INSTR(v.hotc_gl,"Description:")+LENGTH("Description:"))) AS Description,
        SUBSTR(v.hotc_gl,INSTR(v.hotc_gl,"MIM:")+LENGTH("MIM:"),INSTR(v.hotc_gl,"|OMIM_Phenotypes:")-(INSTR(v.hotc_gl,"MIM:")+LENGTH("MIM:"))) AS MIM,
        SUBSTR(v.hotc_gl,INSTR(v.hotc_gl,"OMIM_Phenotypes:")+LENGTH("OMIM_Phenotypes:"),INSTR(v.hotc_gl,"|Orpha_number:")-(INSTR(v.hotc_gl,"OMIM_Phenotypes:")+LENGTH("OMIM_Phenotypes:"))) AS OMIM,
        CASE SUBSTR(v.hotc_gl,INSTR(v.hotc_gl,"MIM:")+LENGTH("MIM:"),INSTR(v.hotc_gl,"|OMIM_Phenotypes:")-(INSTR(v.hotc_gl,"MIM:")+LENGTH("MIM:")))
          WHEN 'NA' THEN ''  
          ELSE '=HYPERLINK("https://www.omim.org/entry/' || SUBSTR(v.hotc_gl,INSTR(v.hotc_gl,"MIM:")+LENGTH("MIM:"),INSTR(v.hotc_gl,"|OMIM_Phenotypes:")-(INSTR(v.hotc_gl,"MIM:")+LENGTH("MIM:"))) || '")' 
        END AS linkOMIM,
        SUBSTR(v.hotc_gl,INSTR(v.hotc_gl,"Orpha_number:")+LENGTH("Orpha_number:"),INSTR(v.hotc_gl,"|Orpha_disorder:")-(INSTR(v.hotc_gl,"Orpha_number:")+LENGTH("Orpha_number:"))) AS ORPHANUM,
        SUBSTR(v.hotc_gl,INSTR(v.hotc_gl,"Orpha_disorder:")+LENGTH("Orpha_disorder:"),INSTR(v.hotc_gl,"|ENSG:")-(INSTR(v.hotc_gl,"Orpha_disorder:")+LENGTH("Orpha_disorder:"))) AS ORPHA,
        v.hotc_panel AS Panel,
        v.transcript AS Ensembl_Transcript,
        v.impact_so AS Consequence,
        v.impact_severity AS Impact_severity,
        REPLACE(v.cdna_position,'/','_') AS Cdna_position,
        REPLACE(v.cds_position,'/','_') AS Cds_position,
        REPLACE(v.codon_change,'/','_') AS Codon_change,
        REPLACE(v.aa_change,'/','_') AS AA_change,

        CAST(r.Variant_impacts_num_refseq AS INT) AS Count_refseq_impacts, 
         r.Grouped_variant_impacts_refseq AS Refseq_impacts, 
        CAST(r.Variant_impacts_num_enst AS INT) AS Count_ensembl_impacts,
         r.Grouped_variant_impacts_enst AS Ensembl_impacts,

        CASE v.is_canonical
          WHEN 1 THEN 'yes'
          WHEN 0 THEN 'no'
        END AS Canonical,
        CASE v.is_splicing
          WHEN 1 THEN 'yes'
          WHEN 0 THEN 'no'
        END AS Splicing,
        CASE v.is_lof
          WHEN 1 THEN 'yes'
          WHEN 0 THEN 'no'
        END AS Loss_of_function,
        v.hc_cadd_phred AS Cadd_phred,
        v.sift_pred AS Sift_pred,
        v.sift_score AS Sift_score,
        v.polyphen_pred AS Polyphen_pred,
        v.polyphen_score AS Polyphen_score,
	CASE
          WHEN LENGTH(v.swissprot)>1 THEN '=HYPERLINK("https://www.uniprot.org/uniprot/' || v.swissprot || '")' 
          ELSE ''
        END AS Swissprot,
        v.hgvsc AS Hgvsc,
        v.hgvsp AS Hgvsp,
	v.hotc_rs_ids AS dbSNP_ID,
        v.hotc_total_ac AS Gnomad_total_ac,
        v.hotc_total_an AS Gnomad_total_an,
        v.hotc_total_af AS Gnomad_total_af,
        v.hotc_ge_ac AS Gnomad_exome_ac,
        v.hotc_ge_an AS Gnomad_exome_an,
        v.hotc_ge_af AS Gnomad_exome_af,
        v.hotc_ge_nhomalt AS Gnomad_exome_nhomalt,
        v.hotc_ge_ac_popmax AS Gnomad_exome_popmax_ac,
        v.hotc_ge_an_popmax AS Gnomad_exome_popmax_an,
        v.hotc_ge_af_popmax AS Gnomad_exome_popmax_af,
        v.hotc_ge_popmax AS Gnomad_exome_popmax,
        v.hotc_ge_af_amr AS Gnomad_exome_af_amr,
        v.hotc_ge_af_fin AS Gnomad_exome_af_fin,
        v.hotc_ge_af_nfe AS Gnomad_exome_af_nfe,
        v.hotc_gg_ac AS Gnomad_genome_ac,
        v.hotc_gg_an AS Gnomad_genome_an,
        v.hotc_gg_af AS Gnomad_genome_af,
        v.hotc_gg_nhomalt AS Gnomad_genome_nhomalt,
        v.hotc_gg_af_amr AS Gnomad_genome_af_amr,
        v.hotc_gg_af_fin AS Gnomad_genome_af_fin,
        v.hotc_gg_af_nfe AS Gnomad_genome_af_nfe,
        CASE 
          WHEN v.hotc_ge_af>0 THEN '=HYPERLINK("https://gnomad.broadinstitute.org/variant/' || REPLACE(v.chrom,'chr','') || '-' || CAST(v.start+1 AS TEXT) || '-' || v.ref || '-' || v.alt || '?dataset=gnomad_r2_1' || '")' 
          WHEN (v.hotc_ge_af <0) AND (v.hotc_gg_af >0) THEN '=HYPERLINK("https://gnomad.broadinstitute.org/variant/' || REPLACE(v.chrom,'chr','') || '-' || CAST(v.start+1 AS TEXT) || '-' || v.ref || '-' || v.alt || '?dataset=gnomad_r3'|| '")' 
          ELSE ''
        END AS linkGNOMAD,
        v.hc_1000gp3_ac AS G1000p3_ac,
        v.hc_1000gp3_af AS G1000p3_af,
        v.hc_1000gp3_afr_ac AS G1000p3_afr_ac,
        v.hc_1000gp3_afr_af AS G1000p3_afr_af,
        v.hc_1000gp3_amr_ac AS G1000p3_amr_ac,
        v.hc_1000gp3_amr_af AS G1000p3_amr_af,
        v.hc_1000gp3_eas_ac AS G1000p3_eas_ac,
        v.hc_1000gp3_eas_af AS G1000p3_eas_af,
        v.hc_1000gp3_eur_ac AS G1000p3_eur_ac,
        v.hc_1000gp3_eur_af AS G1000p3_eur_af,
        v.hc_1000gp3_sas_ac AS G1000p3_sas_ac,
        v.hc_1000gp3_sas_af AS G1000p3_sas_af,

        v.hc_hg19_chr AS Hg19_chr,
        v.hc_hg19_pos AS Hg19_pos,
        v.hotc_clinvar_disease_name AS Clinvar_disease_name,
        v.hotc_clinvar_sig AS Clinvar_sig,
        v.hc_ancestral_allele AS Ancestral_allele,
        v.hc_gerp_nr AS Gerp_nr,
        v.hc_gerp_rs AS Gerp_rs,
        v.hc_gerp_rs_rankscore AS Gerp_rs_rankscore,
        v.hc_integrated_fitcons_score AS Integrated_fitcons_score,
        v.hc_dann_rankscore AS Dann_rankscore,
        v.hc_dann_score AS Dann_score,
        v.hc_eigen_pc_phred AS Eigen_pc_phred,
        v.hc_eigen_pc_raw AS Eigen_pc_raw,
        v.hc_eigen_pc_raw_rankscore AS Eigen_pc_raw_rankscore,
        v.hc_eigen_phred AS Eigen_phred,
        v.hc_eigen_raw AS Eigen_raw,
        v.hc_eigen_coding_or_noncoding AS Eigen_coding_or_noncoding,
        v.hc_mutationassessor_uniprotid AS Mutationassessor_uniprotid,
        v.hc_mutationassessor_pred AS Mutationassessor_pred,
        v.hc_mutationassessor_score AS Mutationassessor_score,
        v.hc_mutationassessor_score_rankscore AS Mutationassessor_score_rankscore,
        v.hc_mutationassessor_variant AS Mutationassessor_variant,
        v.hc_mutationtaster_aae AS Mutationtaster_aae,
        v.hc_mutationtaster_converted_rankscore AS Mutationtaster_converted_rankscore,
        v.hc_mutationtaster_model AS Mutationtaster_model,
        v.hc_mutationtaster_pred AS Mutationtaster_pred,
        v.hc_mutationtaster_score AS Mutationtaster_score,
        v.hc_revel_rankscore AS Revel_rankscore,
        v.hc_revel_score AS Revel_score,
        v.hc_vest3_rankscore AS Vest3_rankscore,
        v.hc_vest3_score AS Vest3_score,
        v.hc_fathmm_mkl_coding_group AS Fathmm_mkl_coding_group,
        v.hc_fathmm_mkl_coding_pred AS Fathmm_mkl_coding_pred,
        v.hc_fathmm_mkl_coding_rankscore AS Fathmm_mkl_coding_rankscore,
        v.hc_fathmm_mkl_coding_score AS Fathmm_mkl_coding_score,
        v.existing_variation AS Existing_variation,
        v.variant_class AS Variant_class,
        v.ccds AS Ccds,
        v.hc_interpro_domain AS Interpro_domain,
        v.hc_gtex_v6p_gene AS Gtex_v6p_gene,
        v.hc_gtex_v6p_tissue AS Gtex_v6p_tissue
    FROM variants v, grouped_impacts_gts r WHERE 
             v.variant_id=r.variant_id;


#select "Variant fields ",count(*) from variant_fields;
#select "";
DROP TABLE IF EXISTS grouped_impacts_refseq;
DROP TABLE IF EXISTS grouped_impacts_enst;
DROP TABLE IF EXISTS var_ids;
DROP TABLE IF EXISTS grouped_impacts_r;
DROP TABLE IF EXISTS grouped_impacts;
DROP TABLE IF EXISTS grouped_impacts_gts;
.mode tabs
.headers on
SELECT * FROM variant_fields;

#DROP TABLE IF EXISTS variant_fields;
# TO TRY LATER 'http://www.mutationtaster.org/cgi-bin/MutationTaster/MutationTaster69.cgi?transcript_stable_id_text=' || v.transcript || '&position_be=' || CAST(v.start+1 AS TEXT) || '&gene=' || v.ensembl_gene_id || '&transcript_stable_id_radio=' || v.transcript || '&sequence_type=CDS&new_base=' || v.alt AS linkMUTST,
