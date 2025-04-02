process COMPUTE_QC {
    input:
    tuple val(meta),
          path(align_flagstat),
          path(markdup_flagstat),
          path(filtered_flagstat),
          path(peak_file),
          path(intersect_file)

    output:
    path "qc_metrics.tsv"       , emit: qc

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python

    def num_of_frags(fl):
        with open(fl, 'r') as f:
            for line in f:
                if "read1" in line:
                    read1_count = int(line.split("+")[0].strip())
                    return read1_count

    def mapping_rate(fl):
        with open(fl, 'r') as f:
            for line in f:
                if "mapped (" in line:
                    start = line.find('(') + 1
                    end = line.find('%')
                    mapping_rate = float(line[start:end]) / 100
                    return mapping_rate

    def dup_rate(fl):
        with open(fl, "r") as file:
            lines = file.readlines()
        duplicates_line = lines[4] 
        mapped_line = lines[6]      

        duplicates = int(duplicates_line.split("+")[0].strip())
        mapped = int(mapped_line.split("+")[0].strip())

        duplicate_rate = duplicates / mapped
        return duplicate_rate

    def count_peak_number(macs2_peaks_file):
        with open(macs2_peaks_file, 'r') as file:
            return sum(1 for _ in file)

    def compute_coverage(sicer_peaks_file):
        total_length = 0
        with open(sicer_peaks_file, 'r') as file:
            for line in file:
                columns = line.strip().split()
                start, end = int(columns[1]), int(columns[2]) 
                total_length += end - start  

        coverage = (total_length * 100) / (3*10**9)  
        coverage = str(round(coverage, 2))+'%'  
        return coverage
    
    def count_reads_in_peaks(bedtools_output_file):
        total_reads = 0
        with open(bedtools_output_file, 'r') as file:
            for line in file:
                total_reads += int(line.strip().split()[-1])  
        return total_reads

    if "$params.broad_peak" != "true":
        header = [
            'sample_id',
            'num_sequenced_fragments',
            'fraction_mappable',
            'num_unique_fragments',
            'non-redundant_fraction',
            'peak_number',
            'frip'
        ]

        content = [
            '${meta.id}',
            str(num_of_frags("${align_flagstat}")),
            str(mapping_rate("${align_flagstat}")),
            str(num_of_frags("${filtered_flagstat}")),
            str(1 - dup_rate("${markdup_flagstat}")),
            str(count_peak_number("${peak_file}")),
            str((count_reads_in_peaks("${intersect_file}")/2/num_of_frags("${filtered_flagstat}")))
        ]
    else:
        header = [
            'sample_id',
            'num_sequenced_fragments',
            'fraction_mappable',
            'num_unique_fragments',
            'non-redundant_fraction',
            'peak_coverage',
            'frip'
        ]

        content = [
            '${meta.id}',
            str(num_of_frags("${align_flagstat}")),
            str(mapping_rate("${align_flagstat}")),
            str(num_of_frags("${filtered_flagstat}")),
            str(1 - dup_rate("${markdup_flagstat}")),
            str(compute_coverage("${peak_file}")),
            str((count_reads_in_peaks("${intersect_file}")/2/num_of_frags("${filtered_flagstat}")))
        ]

    with open('qc_metrics.tsv', 'w') as f:
        f.write('\\t'.join(header) + '\\n')
        f.write('\\t'.join(content))
    """
}

process MERGE_QC {
    input:
    path(files, stageAs: "?/*")

    output:
    path "qc_metrics.tsv"       , emit: qc

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    files = ${files.collect({ "\"${it}\"" })}

    qc = pd.concat([pd.read_csv(f, sep='\\t') for f in files], ignore_index=True)
    qc.to_csv('qc_metrics.tsv', sep='\\t', index=False)
    """
}

workflow REPORT_QC {
    take:
    qc_input

    main:
    COMPUTE_QC ( qc_input ) | collect | MERGE_QC
}