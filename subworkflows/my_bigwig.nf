process MY_GENOMECOV {
    label "process_high"
    container "community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4"
    publishDir "${params.outdir}/bedgraph", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(flagstat)

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    tuple val(meta), path("*.txt")     , emit: scale_factor

    script:
    def prefix = "${meta.id}"
    """
    SCALE_FACTOR=\$(grep '[0-9] mapped (' $flagstat | awk '{print 10000000/\$1}')
    echo \$SCALE_FACTOR > ${prefix}.scale_factor.txt

    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -scale \$SCALE_FACTOR \\
        -pc > tmp.bg

    bedtools sort -i tmp.bg > ${prefix}.bedGraph
    """
}

process MY_BEDGRAPH2BW {
    label "process_high"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:445--h954228d_0' :
        'biocontainers/ucsc-bedgraphtobigwig:445--h954228d_0' }"
    publishDir "${params.outdir}/report/bigwig", mode: 'copy'

    input:
    tuple val(meta), path(bedgraph)
    path sizes

    output:
    tuple val(meta), path("*.bigWig"), emit: bigwig

    script:
    def prefix = "${meta.id}"
    """
    bedGraphToBigWig \\
        $bedgraph \\
        $sizes \\
        ${prefix}.bigWig
    """

}

workflow MY_BIGWIG {
    take:
    ch_bigwig_input
    ch_chrom_sizes

    main:
    MY_GENOMECOV(ch_bigwig_input)
    ch_bedgraph = MY_GENOMECOV.out.bedgraph
    MY_BEDGRAPH2BW(ch_bedgraph,ch_chrom_sizes)
}