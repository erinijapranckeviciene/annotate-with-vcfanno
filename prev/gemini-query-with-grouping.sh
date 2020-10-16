#!/bin/bash
#  gemini query basis from https://github.com/naumenko-sa/cre/blob/master/cre.sh lines 154, 155
#  Sergey's comment: 
#  v.depth = 'None' see https://github.com/chapmanb/bcbio-nextgen/issues/1894

################################################################################################
#  Modification made here : self sufficient query to  extract all required columns
################################################################################################

# Input file
file=$1

sQuery1="select \
        v.chrom as Chrom_geminiv,\
        v.start+1 as Pos_geminiv,\
        v.ref as Ref_geminiv,\
        v.alt as Alt_geminiv,\
        v.impact as Variation_geminiv,\
        v.dp as Depth_dp_geminiv,\
        v.nr as Depth_nr_geminiv,\
        v.nf as Depth_nf_geminiv,\
        v.tc as Depth_tc_geminiv,\
        v.tcf as Depth_tcf_geminiv,\
        v.tcr as Depth_tcr_geminiv,
        v.tr  as Depth_tr_geminiv,\
        v.qual as Quality_geminiv,\
        v.gene as Gene_geminiv,\
        v.clinvar_pathogenic as Clinvar_geminiv,\
        v.cheo_entrezid as Entrezid_cheo,\
        v.cheo_geneinfo as Geneinfo_cheo,\
        v.cheo_mim as Mim_cheo,\
        v.cheo_omimp as Omim_cheo,\
        v.cheo_orpha as Orpha_cheo,\
        v.cheo_myopathygene as Myopathyg_cheo,\
        v.cheo_exon as Exon_names_cheo,\
        v.transcript as Ensembl_transcript_id_geminiv,\
        v.aa_length as AA_position_geminiv,\
        v.exon as Exon_geminiv,\
        v.domains as Protein_domains_geminiv,\
        v.rs_ids as rsIDs_geminiv,\
        v.cheo_af_total as Gnomad_af_cheo,\
        v.cheo_af_popmax as Gnomad_af_popmax_cheo,\
        v.cheo_ac_total as Gnomad_ac_cheo,\
        v.cheo_an_total as Gnomad_an_cheo,\
        v.sift_score as Sift_score_geminiv,\
        v.polyphen_score as Polyphen_score_geminiv,\
        v.cheo_cadd_phred as Cadd_score_cheo,\
        v.cheo_vest3_score as Vest3_score_cheo,\
        v.cheo_revel_score as Revel_score_cheo,\
        v.cheo_gerp_score as Gerp_score_cheo,\
        v.aa_change as AA_change_geminiv,\
        v.hgvsc as Codon_change_geminiv,\
        v.cheo_gerp_score as Gerp_score_cheo,\
        v.cheo_phastcons20way_mammalian as Phastcons20way_mammalian_cheo,\
        v.cheo_phylop20way_mammalian as Conserved_in_20_mammals_cheo,\
        v.hgvsc as Nucleotide_change_ensembl_geminiv,\
	v.hgvsp as Protein_change_ensembl_geminiv,\
        v.old_multiallelic as Old_multiallelic_geminiv,\
        v.callers as Callers_geminiv,\
        (gts).(*),\
        (gt_alt_depths).(*),\
        (gt_depths).(*),\
        (gt_types).(*),\
        count(*) as number_of_impacts,\
        group_concat(i.gene,'|') as gene_geminii,\
        group_concat(i.transcript,'|') as transcript_geminii,\
        group_concat(i.is_exonic,'|') as is_exonic_geminii,\
        group_concat(i.is_coding,'|') as is_coding_geminii,\
        group_concat(i.exon,'|') as exon_geminii,\
        group_concat(i.codon_change,'|') as codon_change_geminii,\
        group_concat(i.aa_change,'|') as aa_change_geminii,\
        group_concat(i.aa_length,'|') as aa_length_geminii,\
        group_concat(i.biotype,'|') as biotype_geminii,\
        group_concat(i.impact,'|') as impact_geminii,\
        group_concat(i.impact_so,'|') as impact_so_geminii,\
        group_concat(i.impact_severity,'|') as impact_severity_geminii,\
        group_concat(i.ccds,'|') as ccds_geminii,\
        group_concat(i.hgvsc,'|') as hgvsc_geminii,\
        group_concat(i.hgvsp,'|') as hgvsp_geminii,\
        group_concat(i.maxentscan_alt,'|') as maxentscan_alt_geminii,\
        group_concat(i.maxentscan_diff,'|') as maxentscan_diff_geminii,\
        group_concat(i.maxentscan_ref,'|') as maxentscan_ref_geminii,\
        group_concat(i.spliceregion,'|') as spliceregion_geminii,\
        group_concat(i.polyphen_pred,'|') as polyphen_pred_geminii,\
        group_concat(i.polyphen_score,'|') as polyphen_score_geminii,\
        group_concat(i.sift_pred,'|') as sift_pred_geminii,\
        group_concat(i.sift_score, '|') as sift_score_geminii,\
        group_concat( i.gene || '; exon: ' || i.exon || ';' || i.hgvsc || ';' || i.hgvsp, '|')  as info_geminii,\
        v.variant_id as Variant_id_v from variants v, variant_impacts i where\
 ( (v.variant_id=i.variant_id) and\
   ( (i.impact_severity in ('HIGH', 'MED')) and (v.impact_severity in ('HIGH', 'MED') ) ) and\
   ( (v.cheo_af_popmax <=0.01) or (v.cheo_af_popmax >=99.9) ) and\
   ( (v.dp >=10) or (v.nr>=10) ) ) group by v.variant_id"


#echo ${sQuery1}
#echo ${file}

    gemini query --header -q "${sQuery1}" $file

