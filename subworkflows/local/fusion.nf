include { ARRIBA_ARRIBA                       } from '../../modules/nf-core/arriba/arriba/main'
include { STAR_ALIGN as STAR_FOR_STARFUSION   } from '../../modules/nf-core/star/align/main'
include { STARFUSION                          } from '../../modules/local/starfusion/detect/main'
include { FUSIONCATCHER_DETECT                } from '../../modules/local/fusioncatcher/detect/main'
include { ONCOKB_FUSIONANNOTATOR              } from '../../modules/local/oncokb/fusionannotator/main'
include { AGFUSION_BATCH                      } from '../../modules/local/agfusion/batch/main'
include { AGFUSION_BATCH as AGFUSION_CLINICAL } from '../../modules/local/agfusion/batch/main'
include { TO_CFF as ARRIBA_TO_CFF             } from '../../modules/local/convert_to_cff/main'
include { TO_CFF as FUSIONCATCHER_TO_CFF      } from '../../modules/local/convert_to_cff/main'
include { TO_CFF as STARFUSION_TO_CFF         } from '../../modules/local/convert_to_cff/main'
include { CAT_CAT as MERGE_CFF                } from '../../modules/nf-core/cat/cat/main'
include { METAFUSION_RUN                      } from '../../modules/local/metafusion/run/main'
include { ADD_FLAG                            } from '../../modules/local/add_flags/main'
include { CFF_ANNOTATE as CFF_FINALIZE        } from '../../modules/local/cff_annotate/main'
include { CFF_ANNOTATE as ADD_FLAG_AGFUSION   } from  '../../modules/local/cff_annotate/main'
include { FUSION_FILTER                       } from  '../../modules/local/fusion_filtering/main'
include { ANNOTATION_RUN                     } from '../../modules/local/annotation/run/main'
include { REFS_ANNOTATIONS                     } from '../../modules/local/annotation/refs/main'

workflow FUSION {

    take:
    reads
    reads_untrimmed
    bam
    star_index
    fasta
    gtf
    starfusion_ref
    fusioncatcher_ref
    agfusion_db
    pyensembl_cache
    gene_bed
    gene_info
    blocklist
    arriba_blacklist
    arriba_known_fusions
    arriba_protein_domains
    clinical_genes
    transcript_allowlist

    main:
    ch_versions = Channel.empty()
    //fasta = params.fasta
    //gene_bed = params.metafusion_gene_bed
    //blocklist = params.metafusion_blocklist

    ARRIBA_ARRIBA(
        bam,
        fasta,
        gtf,
        arriba_blacklist.map{[[:],it]},
        arriba_known_fusions.map{[[:],it]},
        [[:],[]],
        [[:],[]],
        arriba_protein_domains.map{[[:],it]}
    )
    ch_versions = ch_versions.mix(ARRIBA_ARRIBA.out.versions.first())

    STAR_FOR_STARFUSION(
        reads,
        // use the star index in the starfusion reference to ensure compatibility
        starfusion_ref.map{ [[id:params.genome],file( it + "/ref_genome.fa.star.idx")] },
        starfusion_ref.map{ [[id:params.genome],file( it + "/ref_annot.gtf")] },
        false,
        [],
        []
    )
    ch_versions = ch_versions.mix(STAR_FOR_STARFUSION.out.versions.first())

    //reads_junction = reads.join( STAR_FOR_STARFUSION.out.junction,by:[0] )

    STARFUSION(
        STAR_FOR_STARFUSION.out.junction.map{ meta, junction -> [ meta, [], junction ] },
        starfusion_ref
    )
    ch_versions = ch_versions.mix(STARFUSION.out.versions.first())

    FUSIONCATCHER_DETECT(
        reads_untrimmed,
        fusioncatcher_ref
    )
    ch_versions = ch_versions.mix(FUSIONCATCHER_DETECT.out.versions.first())

    fc_fusions = ["GRCh37","hg19","smallGRCh37"].contains(params.genome) ? FUSIONCATCHER_DETECT.out.fusions_alt : FUSIONCATCHER_DETECT.out.fusions


    ARRIBA_TO_CFF(ARRIBA_ARRIBA.out.fusions
            .map{ meta, file ->[ meta, "arriba", file ] })
    FUSIONCATCHER_TO_CFF(fc_fusions
                    .map{ meta, file -> [ meta, "fusioncatcher", file ] } )
    STARFUSION_TO_CFF(STARFUSION.out.abridged
                    .map{ meta, file -> [ meta, "starfusion", file ] })
    // get expected number of callers for groupTuple
    numcallers = 1 + ( params.starfusion_url ? 1 : 0 ) + ( ["GRCh37","GRCh38"].contains(params.genome) ? 1 : 0 )

    MERGE_CFF(
        ARRIBA_TO_CFF.out.cff
            .map{ meta, file -> [meta, file]}
            .mix(
                FUSIONCATCHER_TO_CFF.out.cff
                    .map{ meta, file -> [meta, file]}
            ).mix(
                STARFUSION_TO_CFF.out.cff
                    .map{ meta, file -> [meta, file]}
            ).groupTuple(by:[0],size:numcallers),
    )

    METAFUSION_RUN(
        MERGE_CFF.out.file_out,
        gene_bed.map{ it[1] }.first(),
        gene_info.map{ it[1] }.first(),
        fasta.map{ it[1] }.first(),
        blocklist,
        transcript_allowlist
    )

    ADD_FLAG(
        METAFUSION_RUN.out.cluster
            .join(METAFUSION_RUN.out.cis, by:0)
            .join(METAFUSION_RUN.out.cff, by:0)
            .join(METAFUSION_RUN.out.problem_chrom, by:0)
            .join(METAFUSION_RUN.out.filters, by:0),
        clinical_genes
    )

    ONCOKB_FUSIONANNOTATOR(ADD_FLAG.out.unfiltered_cff)
    ch_versions = ch_versions.mix(ONCOKB_FUSIONANNOTATOR.out.versions.first())

    AGFUSION_BATCH(
        ADD_FLAG.out.unfiltered_cff,
        agfusion_db,
        pyensembl_cache
    )
    ch_versions = ch_versions.mix(AGFUSION_BATCH.out.versions.first())

    if (params.run_oncokb_fusionannotator) {
        CFF_FINALIZE(
            ADD_FLAG.out.unfiltered_cff
                .join(ONCOKB_FUSIONANNOTATOR.out.oncokb_fusions, by:0)
                .join(AGFUSION_BATCH.out.fusion_transcripts_tsv, by:0)
                .map{ meta, cff, oncokb, agfusion_file ->
                    [ meta, cff, oncokb, agfusion_file ]
                },
            transcript_allowlist
        )
    } else {
        CFF_FINALIZE(
            ADD_FLAG.out.unfiltered_cff
                .join(AGFUSION_BATCH.out.fusion_transcripts_tsv, by:0)
                .map{ meta, cff, agfusion_file ->
                    [ meta, cff, [], agfusion_file ]
                },
            transcript_allowlist
        )
    }

    FUSION_FILTER(
        CFF_FINALIZE.out.filtered_cff
             .join(STARFUSION.out.coding_effect, by:0)
             .join(FUSIONCATCHER_DETECT.out.fusions, by:0)
             .join(ARRIBA_ARRIBA.out.fusions, by:0),
        clinical_genes
    )

    REFS_ANNOTATIONS(
        params.genome,
        params.ensembl_version,
        params.baits.idt_v2.baits,
        params.gtf,
        params.genomes.(params.genome).ucsc_Pfam
    )

    ANNOTATION_RUN(
        FUSION_FILTER.out.filtered_fusions,
        REFS_ANNOTATIONS.out.reference_transcripts,
        REFS_ANNOTATIONS.out.sv_table,
        REFS_ANNOTATIONS.out.refFlat,
        REFS_ANNOTATIONS.out.refFlat_summary,
        REFS_ANNOTATIONS.out.kinase_domains,
        params.reference_base + '/annotation/tumourSuppressors_IMPACT.txt',
        params.reference_base + '/annotation/oncokb_known_fusions.txt',
        params.reference_base + '/annotation/keygenes.txt'
    )

    AGFUSION_CLINICAL(
        ADD_FLAG.out.unfiltered_clinical_cff,
        agfusion_db,
        pyensembl_cache
    )

    ADD_FLAG_AGFUSION(
        ADD_FLAG.out.unfiltered_cff
            .join(AGFUSION_CLINICAL.out.fusion_transcripts_tsv, by:0)
            .map{ meta, cff, agfusion_file ->
                    [ meta, cff, [], agfusion_file]
                },
        transcript_allowlist
    )

    ch_versions = ch_versions.mix(CFF_FINALIZE.out.versions.first())
    ch_versions = ch_versions.mix(FUSION_FILTER.out.versions.first())
    ch_versions = ch_versions.mix(AGFUSION_CLINICAL.out.versions.first())
    ch_versions = ch_versions.mix(ADD_FLAG_AGFUSION.out.versions.first())
    ch_versions = ch_versions.mix(ADD_FLAG.out.versions.first())
    ch_versions = ch_versions.mix(METAFUSION_RUN.out.versions.first())
    ch_versions = ch_versions.mix(ARRIBA_TO_CFF.out.versions.first())
    ch_versions = ch_versions.mix(FUSIONCATCHER_TO_CFF.out.versions.first())
    ch_versions = ch_versions.mix(STARFUSION_TO_CFF.out.versions.first())
    ch_versions = ch_versions.mix(ADD_FLAG_AGFUSION.out.versions.first())


    emit:
    ch_versions
}
