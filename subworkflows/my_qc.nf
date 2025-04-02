process MY_DEDUP {
    label "process_high"
    container "community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c"
    publishDir "${params.outdir}/dedup", mode: 'copy'

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("*.target.deduplicated.sorted.bam"), emit: bam
    tuple val(meta), path("*.flagstat"), emit: flagstat

    script:
    """
    bash_bam=${bam_file}
    prefix="\${bash_bam/.target.markdup.sorted.bam/}"
    dedup_bam_file="\${prefix}.target.deduplicated.sorted.bam"
    flagstat_file="\${prefix}_filtered.flagstat"
    samtools view -@ 32 -b -F 1024 "$bam_file" > "\$dedup_bam_file"
    samtools flagstat "\$dedup_bam_file" > "\$flagstat_file"
    """

}


process BAM2BED {
    label "process_high"
    container "community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4"
    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("*.bed"), emit: bed

    script:
    """
    bash_bam=${bam_file}
    prefix="\${bash_bam/.target.deduplicated.sorted.bam/}"
    bedtools bamtobed -i "\$bash_bam" > "\${prefix}.bed"
    """
}

process CALL_NARROW_PEAK{
    label "process_high"
    container "community.wave.seqera.io/library/macs3:3.0.3--22198510e77bb7aa"
    publishDir "${params.outdir}/narrow_peak", mode: 'copy'

    input:
    tuple val(meta), path(bed_file)

    output:
    tuple val(meta), path("*.narrowPeak"), emit: narrowPeak

    script:
    """
    bash_bed=${bed_file}
    prefix="\${bash_bed/.bed/}"
    macs3 callpeak -t "$bed_file" -f BED -g hs -n "\$prefix" --nomodel --shift -100 --extsize 200 -q 0.05  --keep-dup all
    """
}

process CALL_BROAD_PEAK{
    label "process_high"
    container "community.wave.seqera.io/library/sicer2:1.0.3--0c8811a485e8a75e"
    publishDir "${params.outdir}/broad_peak", mode: 'copy'

    input:
    tuple val(meta), path(bed_file)

    output:
    tuple val(meta), path("*.scoreisland"), emit: scoreisland

    script:
    """
    sicer -w 1000 -rt 16 -f 300 -egf 0.8 -fdr 0.01 -g 3000 -e 100 -s hg38 -t "$bed_file" 
    """
}

process INTERSECT {
    label "process_high"
    publishDir "${params.outdir}/bedtools_intersect", mode: 'copy'
    container "community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4"

    input:
    tuple val(meta), path(bam_file), path(peak_file)

    output:
    tuple val(meta), path("*_bedtools.txt"), emit: intersect

    script:
    """
    bash_bam=${bam_file}
    prefix="\${bash_bam/.target.deduplicated.sorted.bam/}"
    bedtools intersect -a "$bam_file" -b "$peak_file" -c -bed > "\${prefix}_bedtools.txt"
    """
}

workflow MY_QC {
    take:
    ch_samtools_bam

    main:
    MY_DEDUP(ch_samtools_bam)
    ch_deduped_bam = MY_DEDUP.out.bam

    BAM2BED(ch_deduped_bam)
    ch_bed = BAM2BED.out.bed

    if (params.broad_peak == false){
        CALL_NARROW_PEAK(ch_bed)
        ch_intersect_peak = CALL_NARROW_PEAK.out.narrowPeak
        ch_output_peak = CALL_NARROW_PEAK.out.narrowPeak
    } else {
        CALL_BROAD_PEAK(ch_bed)
        ch_intersect_peak = CALL_BROAD_PEAK.out.scoreisland
        ch_output_peak = CALL_BROAD_PEAK.out.scoreisland
    }
    ch_intersect_bam = MY_DEDUP.out.bam

    ch_intersect_bam
    .combine(ch_intersect_peak, by: 0)
    .set { ch_intersect_input }
    INTERSECT(ch_intersect_input)

    emit:
    flagstat = MY_DEDUP.out.flagstat
    peak = ch_output_peak
    intersect = INTERSECT.out.intersect

}