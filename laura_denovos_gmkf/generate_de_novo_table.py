from collections import defaultdict, Counter
from gnomad_hail import *
from gnomad_hail.resources.sample_qc import *
from gnomad_hail.utils.generic import file_exists, filter_low_conf_regions
import hail as hl
import os
import subprocess

print(subprocess.check_output(
    "/opt/conda/default/bin/python3 -m pip install git+https://github.com/bw2/hail-utils.git", 
    shell=True, encoding="UTF-8"))

from hail_utils.vep import compute_derived_vep_fields, get_expr_for_vep_sorted_transcript_consequences_array, get_expr_for_vep_consequence_terms
from hail_utils.info_fields import recompute_AC_AN_AF, annotate_in_LCR, annotate_in_segdup
from hail_utils.filters import filter_out_variants_where_all_samples_are_hom_ref, filter_out_LCRs, filter_out_segdups, filter_to_autosomes
from hail_utils.variants import get_expr_for_variant_type
from hail_utils.io import file_exists
from hail_utils.infer_pedigree import impute_sex, compute_kinship_ht, infer_families, get_duplicated_samples_ibd
from hail_utils.denovos import compute_mendel_denovos, compute_samocha_denovos

import argparse

p = argparse.ArgumentParser()
p.add_argument("--force", action="store_true")
p.add_argument("--mendelian", action="store_true")
p.add_argument("data_label", nargs="+", choices=["pcr_free", "pcr_plus", "rgp"])
args = p.parse_args()

BASE_DIR = "gs://seqr-bw2/project__laura_denovo_vcfs__gmkf"
PCR_FREE_VCF_PATH = os.path.join(BASE_DIR, "GMKF_Gabriel_Luquetti_Craniofacial_WGS_v2.vcf.gz")
PCR_PLUS_VCF_PATH = os.path.join(BASE_DIR, "GMKF_Gabriel_Luquetti_Craniofacial_WGS_PCR_Plus_v2.vcf.gz")
RGP_VCF_PATH = "gs://seqr-datasets/GRCh38/RDG_WGS_Broad_Internal/v7/sharded_vcf/RDG_WGS_Broad_Internal.filtered.*.vcf.gz"

force = args.force
algo_label = "mendel" if args.mendelian else "samocha"

for data_label in args.data_label:
    print(f"Processing {data_label}. Force: {force}")
    if data_label == "pcr_plus":
        VCF_PATH = PCR_PLUS_VCF_PATH
    elif data_label == "pcr_free":
        VCF_PATH = PCR_FREE_VCF_PATH
    elif data_label == "rgp":
        VCF_PATH = RGP_VCF_PATH
    else:
        raise ValueError("Invalid data label: " + data_label)

    mt_path = os.path.join(BASE_DIR, f"{data_label}.mt")
    print(mt_path)

    filtered_mt = os.path.join(BASE_DIR, f"{data_label}.basic_filters.mt") 
    if force or not file_exists(filtered_mt):
        print(f"{data_label}: Generating {filtered_mt}")
        mt = hl.import_vcf(VCF_PATH, force_bgz=True, min_partitions=10000, reference_genome="GRCh38")

        if data_label == "rgp":
            samples_to_keep = hl.set({s for s in mt.s.collect() if "RGP" in s.upper()})
            mt = mt.filter_cols(samples_to_keep.contains(mt['s']), keep=True)

        mt = filter_out_variants_where_all_samples_are_hom_ref(mt)
        mt = hl.split_multi_hts(mt)
        mt = recompute_AC_AN_AF(mt)
        mt = filter_to_adj(mt)

        #mt = filter_out_LCRs(mt)
        #mt = filter_out_segdups(mt)
        mt = annotate_in_LCR(mt)
        mt = annotate_in_segdup(mt)
    
        mt = mt.checkpoint(filtered_mt, overwrite=force, _read_if_exists=not force)

    else:
        print(f"{data_label}: Reading {filtered_mt}...")
        mt = hl.read_matrix_table(filtered_mt)

    total_variants, num_samples = mt.count()
    print(f"{num_samples} VCF samples")
    print(f"{total_variants} variants")

    print(f"{data_label}: impute sex...")
    vcf_samples = mt.s.collect()
    imputed_sex = impute_sex(mt)

    file_path = os.path.join(BASE_DIR, f"{data_label}.kinship.ht")  
    if not file_exists(file_path):
        print(f"{data_label}: Generating {file_path}")
        kin_ht = compute_kinship_ht(mt)
        kin_ht = kin_ht.checkpoint(file_path, overwrite=force, _read_if_exists=not force)
    else:
        print(f"Reading table {file_path}")
        kin_ht = hl.read_table(file_path)
    
    print(f"{data_label}: infer families...")
    pedigree, duos, decisions = infer_families(kin_ht, imputed_sex)
        
    pedigree = hl.Pedigree(pedigree.complete_trios())
    pedigree.write(os.path.join(BASE_DIR, f"{data_label}.pediree.fam"))
    pedigree.complete_trios()

    print(len(pedigree.complete_trios())*3)
    print(len(pedigree.complete_trios())*3/len(vcf_samples))

    file_path = os.path.join(BASE_DIR, f"{data_label}.{algo_label}_denovos.ht")
    if args.mendelian:
        if not file_exists(file_path):
            print(f"{data_label}: Generating {file_path}")
            denovos = compute_mendel_denovos(mt, pedigree)
            denovos = denovos.checkpoint(file_path, overwrite=True, _read_if_exists=not force)
        else:
            print(f"Reading table {file_path}")
            denovos = hl.read_table(file_path)

    else:
        if not file_exists(file_path):
            print(f"{data_label}: Generating {file_path}")
            denovos = compute_samocha_denovos(mt, pedigree)
            denovos = denovos.checkpoint(file_path, overwrite=True, _read_if_exists=not force)
        else:
            print(f"Reading table {file_path}")
            denovos = hl.read_table(file_path)

    file_path = os.path.join(BASE_DIR, f"{data_label}.{algo_label}_denovos.vep.ht")
    if not file_exists(file_path):
        print(f"{data_label}: Generating {file_path}")
        denovos = denovos.key_by('locus', 'alleles')
        denovos = hl.vep(denovos, "file:///vep_data/vep-gcloud.json", name="vep", block_size=100)

        denovos = denovos.checkpoint(file_path, overwrite=True, _read_if_exists=not force)
    else:
        print(f"Reading table {file_path}")
        denovos = hl.read_table(file_path)


    print(f"{data_label}: annotate and export {file_path}")
    denovos = denovos.annotate(sorted_transcript_consequences = get_expr_for_vep_sorted_transcript_consequences_array(denovos.vep))
    denovos = denovos.annotate(transcript_consequence_terms = get_expr_for_vep_consequence_terms(denovos.sortedTranscriptConsequences))
    denovos = denovos.annotate(transcript_consequence_categories = denovos.sorted_transcript_consequences.map(lambda c: c.category))
    denovos = denovos.annotate(transcript_consequence = hl.cond(hl.len(denovos.transcript_consequence_terms) > 0, denovos.transcript_consequence_terms[0], "other"))
    denovos = denovos.annotate(transcript_consequence_category = hl.cond(hl.len(denovos.transcript_consequence_categories) > 0, denovos.transcript_consequence_categories[0], "other"))

    joined_mt_rows = mt.rows()[(denovos.locus, denovos.alleles)]

    gnomad_genomes = hl.read_table("gs://seqr-reference-data/GRCh38/gnomad/gnomad.genomes.r2.1.1.sites.liftover_grch38.ht")
    gnomad_exomes = hl.read_table("gs://seqr-reference-data/GRCh38/gnomad/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht")

    joined_gnomad_genomes_ht = gnomad_genomes[(denovos.locus, denovos.alleles)]
    joined_gnomad_exomes_ht = gnomad_exomes[(denovos.locus, denovos.alleles)]

    file_path = os.path.join(BASE_DIR, f"{data_label}.{algo_label}_de_novos_table.tsv")
    denovos = denovos.annotate(
        dataset=data_label,
        variant_type=get_expr_for_variant_type(denovos),
        in_LCR=joined_mt_rows.info.in_LCR,
        in_segdup=joined_mt_rows.info.in_segdup,
        filters=hl.cond(hl.len(joined_mt_rows.filters) > 0, hl.delimit(joined_mt_rows.filters, ','), 'PASS'),
        AC=joined_mt_rows.info.AC,
        AF=joined_mt_rows.info.AF, 
        QD=joined_mt_rows.info.QD,
        gnomAD_genomes_AF = joined_gnomad_genomes_ht.freq[0].AF,
        gnomAD_genomes_AC = joined_gnomad_genomes_ht.freq[0].AC,
        gnomAD_exomes_AF = joined_gnomad_exomes_ht.freq[0].AF,
        gnomAD_exomes_AC = joined_gnomad_exomes_ht.freq[0].AC,
    )

    denovos = denovos.key_by()
    if args.mendelian:
        denovos = denovos.select(
            'dataset',
            's', 'mendel_code',
            'transcript_consequence', 'transcript_consequence_category',
            'variant_type', 'filters', 'AC', 'AF', 'QD', 'in_LCR', 'in_segdup')
    else:
        denovos = denovos.select(
            'dataset',
            'confidence',
            'proband', 'father', 'mother', 'proband_AB', 'proband_DP', 'proband_GQ', 'mother_AB', 'mother_DP', 'mother_GQ', 'father_AB', 'father_DP', 'father_GQ',
            'transcript_consequence', 'transcript_consequence_category',
            'variant_type', 'filters', 'AC', 'AF', 'QD', 'in_LCR', 'in_segdup')


    denovos.export(file_path)

    print(f"Done with {data_label} {algo_label}!")
