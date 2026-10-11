# =============================================================================
# Mite Virome RNA-seq Pipeline
# Species: Neoseiulus californicus (Ca) and Neoseiulus cucumeris (NC)
# Samples: Ca1-3 and NC1-3, paired-end
# =============================================================================

configfile: "config.yaml"

SAMPLES = {
    "Ca": ["Ca1", "Ca2", "Ca3"],
    "NC": ["NC1", "NC2", "NC3"],
}
ALL_SAMPLES = SAMPLES["Ca"] + SAMPLES["NC"]
SPECIES     = list(SAMPLES.keys())
SAMPLE_TO_SPECIES = {s: sp for sp, samples in SAMPLES.items() for s in samples}

# Tool environment activation paths
ENV = {
    "multiqc":    ". /data/bin/miniconda2/envs/multiqc-v1.9/env_multiqc.sh",
    "trimgalore": ". /data/bin/miniconda2/envs/trimGalore-v0.6.6/env_trimGalore.sh",
    "star":       ". /data/bin/miniconda2/envs/star-v2.7.11/env_star.sh",
    "bowtie2":    ". /data/bin/miniconda2/envs/bowtie2-v2.5.4/env_bowtie2.sh",
    "samtools":   ". /data/bin/miniconda2/envs/samtools-v1.22.1/env_samtools.sh",
    "kraken2":    ". /data/bin/miniconda2/envs/kraken2-v2.1.2/env_kraken2.sh",
    "trinity":    ". /data/bin/miniconda2/envs/trinity-v2.15.1/env_trinity.sh",
    "diamond":    ". /data/bin/miniconda2/envs/diamond-v2.1.9/env_diamond.sh;",
    "salmon_env": ". /data/bin/miniconda2/envs/salmon-v1.9.0/env_salmon.sh",  # check exact path
}


# =============================================================================
# Target rule
# =============================================================================
rule all:
    input:
        # QC raw
        expand("results/QC/raw/{sample}_R1_001_fastqc.html", sample=ALL_SAMPLES),
        expand("results/QC/raw/{sample}_R2_001_fastqc.html", sample=ALL_SAMPLES),
        "results/QC/raw/multiqc_report.html",
        # QC trimmed
        expand("results/QC/trimmed/{sample}_R1_001_val_1_fastqc.html", sample=ALL_SAMPLES),
        expand("results/QC/trimmed/{sample}_R2_001_val_2_fastqc.html", sample=ALL_SAMPLES),
        "results/QC/trimmed/multiqc_report.html",
        # Host mapping
        expand("results/host_mapping/{sample}/Log.final.out", sample=ALL_SAMPLES),
        # Unmapped reads
        expand("results/unmapped/{sample}_unmapped_R1.fq.gz", sample=ALL_SAMPLES),
        expand("results/unmapped/{sample}_unmapped_R2.fq.gz", sample=ALL_SAMPLES),
        # rRNA
        expand("results/rrna/{sample}.flagstat", sample=ALL_SAMPLES),
        # Kraken2
        expand("results/kraken2/{sample}.report", sample=ALL_SAMPLES),
        # Trinity
        expand("results/trinity/{species}_trinity.Trinity.fasta", species=SPECIES),
        expand("results/trinity/{species}_trinity.assembly_stats.txt", species=SPECIES),
        # DIAMOND
        expand("results/diamond/{species}/blastx.tsv", species=SPECIES),
        # Taxonomy
        expand("results/taxonomy/{species}/taxonomy_summary.tsv", species=SPECIES),
        # Summary
        "results/summary/pipeline_summary.tsv",
        # STEP 11 — Salmon quantification (New Addition)
        expand("results/salmon/{species}/merged_tpm.tsv", species=SPECIES),
        # STEP 12 — Final Integrated Summary
        expand("results/summary/{species}_viral_abundance_summary.tsv", species=SPECIES),


# =============================================================================
# STEP 1 — Raw QC
# =============================================================================
rule fastqc_raw:
    input:
        r1 = lambda wc: f"{config['raw_data_dir']}/{wc.sample}_R1_001.fastq.gz",
        r2 = lambda wc: f"{config['raw_data_dir']}/{wc.sample}_R2_001.fastq.gz",
    output:
        html_r1 = "results/QC/raw/{sample}_R1_001_fastqc.html",
        html_r2 = "results/QC/raw/{sample}_R2_001_fastqc.html",
        zip_r1  = "results/QC/raw/{sample}_R1_001_fastqc.zip",
        zip_r2  = "results/QC/raw/{sample}_R2_001_fastqc.zip",
    threads: 2
    resources:
        mem_mb  = 4000,
        runtime = 60,
    log:
        "logs/fastqc_raw/{sample}.log"
    shell:
        "fastqc -t {threads} -o results/QC/raw/ {input.r1} {input.r2} 2> {log}"


rule multiqc_raw:
    input:
        expand("results/QC/raw/{sample}_R1_001_fastqc.zip", sample=ALL_SAMPLES),
        expand("results/QC/raw/{sample}_R2_001_fastqc.zip", sample=ALL_SAMPLES),
    output:
        "results/QC/raw/multiqc_report.html"
    log:
        "logs/multiqc_raw.log"
    shell:
        """
        {ENV[multiqc]}
        multiqc results/QC/raw/ -o results/QC/raw/ --filename multiqc_report.html --force 2> {log}
        """


# =============================================================================
# STEP 2 — Trimming + post-trim QC
# =============================================================================
rule trim_galore:
    input:
        r1 = lambda wc: f"{config['raw_data_dir']}/{wc.sample}_R1_001.fastq.gz",
        r2 = lambda wc: f"{config['raw_data_dir']}/{wc.sample}_R2_001.fastq.gz",
    output:
        r1  = "results/trimmed/{sample}_R1_001_val_1.fq.gz",
        r2  = "results/trimmed/{sample}_R2_001_val_2.fq.gz",
        rep = "results/trimmed/{sample}_R1_001.fastq.gz_trimming_report.txt",
    threads: 4
    resources:
        mem_mb  = 8000,
        runtime = 120,
    params:
        extra = config.get("trim_galore_extra", "--quality 20 --length 50"),
    log:
        "logs/trim_galore/{sample}.log"
    shell:
        """
        {ENV[trimgalore]}
        trim_galore --paired --cores {threads} \
            {params.extra} \
            -o results/trimmed/ \
            {input.r1} {input.r2} 2> {log}
        """

rule fastqc_trimmed:
    input:
        r1 = "results/trimmed/{sample}_R1_001_val_1.fq.gz",
        r2 = "results/trimmed/{sample}_R2_001_val_2.fq.gz",
    output:
        html_r1 = "results/QC/trimmed/{sample}_R1_001_val_1_fastqc.html",
        html_r2 = "results/QC/trimmed/{sample}_R2_001_val_2_fastqc.html",
        zip_r1  = "results/QC/trimmed/{sample}_R1_001_val_1_fastqc.zip",
        zip_r2  = "results/QC/trimmed/{sample}_R2_001_val_2_fastqc.zip",
    threads: 2
    resources:
        mem_mb  = 4000,
        runtime = 60,
    log:
        "logs/fastqc_trimmed/{sample}.log"
    shell:
        "fastqc -t {threads} -o results/QC/trimmed/ {input.r1} {input.r2} 2> {log}"


rule multiqc_trimmed:
    input:
        expand("results/QC/trimmed/{sample}_R1_001_val_1_fastqc.zip", sample=ALL_SAMPLES),
        expand("results/QC/trimmed/{sample}_R2_001_val_2_fastqc.zip", sample=ALL_SAMPLES),
    output:
        "results/QC/trimmed/multiqc_report.html"
    log:
        "logs/multiqc_trimmed.log"
    shell:
        """
        {ENV[multiqc]}
        multiqc results/QC/trimmed/ -o results/QC/trimmed/ --filename multiqc_report.html --force 2> {log}
        """


# =============================================================================
# STEP 3 — Build STAR index per species
# =============================================================================
rule star_genome_index:
    input:
        fasta = lambda wc: config["genomes"][wc.species]["fasta"],
    output:
        sa = "results/star_index/{species}/SA",
    threads: 16
    resources:
        mem_mb  = 40000,
        runtime = 120,
    params:
        outdir          = "results/star_index/{species}/",
        sa_index_nbases = lambda wc: config["genomes"][wc.species].get("star_genomeSAindexNbases", 11),
        gtf_flag        = lambda wc: (
            f"--sjdbGTFfile {config['genomes'][wc.species]['gtf']}"
            if config["genomes"][wc.species].get("gtf") else ""
        ),
    log:
        "logs/star_index/{species}.log"
    shell:
        """
        {ENV[star]}
        mkdir -p {params.outdir}
        STAR --runMode genomeGenerate \
            --genomeDir {params.outdir} \
            --genomeFastaFiles {input.fasta} \
            {params.gtf_flag} \
            --genomeSAindexNbases {params.sa_index_nbases} \
            --runThreadN {threads} 2> {log}
        """


# =============================================================================
# STEP 4 — Map to species-specific genome, export unmapped reads
# =============================================================================
rule star_map:
    input:
        r1 = "results/trimmed/{sample}_R1_001_val_1.fq.gz",
        r2 = "results/trimmed/{sample}_R2_001_val_2.fq.gz",
        sa = lambda wc: f"results/star_index/{SAMPLE_TO_SPECIES[wc.sample]}/SA",
    output:
        bam      = "results/host_mapping/{sample}/Aligned.sortedByCoord.out.bam",
        log_file = "results/host_mapping/{sample}/Log.final.out",
        unmap_r1 = "results/host_mapping/{sample}/Unmapped.out.mate1",
        unmap_r2 = "results/host_mapping/{sample}/Unmapped.out.mate2",
    params:
        index_dir = lambda wc: f"results/star_index/{SAMPLE_TO_SPECIES[wc.sample]}/",
        outprefix = "results/host_mapping/{sample}/",
    threads: 16
    resources:
        mem_mb  = 40000,
        runtime = 180,
    log:
        "logs/star_map/{sample}.log"
    shell:
        """
        {ENV[star]}
        STAR --runThreadN {threads} \
            --genomeDir {params.index_dir} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix {params.outprefix} \
            --outTmpDir {params.outprefix}_tmp \
            --limitBAMsortRAM 10000000000 \
            --outFilterMultimapNmax 10 \
            --outFilterMismatchNmax 10 \
            2> {log}
        """


rule compress_unmapped:
    input:
        r1 = "results/host_mapping/{sample}/Unmapped.out.mate1",
        r2 = "results/host_mapping/{sample}/Unmapped.out.mate2",
    output:
        r1 = "results/unmapped/{sample}_unmapped_R1.fq.gz",
        r2 = "results/unmapped/{sample}_unmapped_R2.fq.gz",
    threads: 4
    resources:
        mem_mb  = 4000,
        runtime = 60,
    shell:
        """
        {ENV[samtools]}
        pigz -p {threads} -c {input.r1} > {output.r1}
        pigz -p {threads} -c {input.r2} > {output.r2}
        """


# =============================================================================
# STEP 5 — rRNA screen (SILVA SSU + LSU)
# =============================================================================
rule bowtie2_rrna_index:
    input:
        config["rrna_fasta"]
    output:
        multiext("results/rrna_index/rrna",
                 ".1.bt2", ".2.bt2", ".3.bt2",
                 ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    threads: 8
    resources:
        mem_mb  = 8000,
        runtime = 60,
    log:
        "logs/bowtie2_rrna_index.log"
    shell:
        """
        {ENV[bowtie2]}
        mkdir -p results/rrna_index/
        bowtie2-build --threads {threads} {input} results/rrna_index/rrna 2> {log}
        """


rule map_rrna:
    input:
        r1    = "results/trimmed/{sample}_R1_001_val_1.fq.gz",
        r2    = "results/trimmed/{sample}_R2_001_val_2.fq.gz",
        index = multiext("results/rrna_index/rrna",
                         ".1.bt2", ".2.bt2", ".3.bt2",
                         ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    output:
        bam      = "results/rrna/{sample}.bam",
        flagstat = "results/rrna/{sample}.flagstat",
    threads: 8
    resources:
        mem_mb  = 32000,
        runtime = 120,
    log:
        "logs/rrna/{sample}.log"
    shell:
        """
        {ENV[bowtie2]}
        {ENV[samtools]}
        bowtie2 -p {threads} --very-sensitive \
            -x results/rrna_index/rrna \
            -1 {input.r1} -2 {input.r2} \
            -S results/rrna/{wildcards.sample}.sam \
            2> {log}
        samtools sort -@ {threads} \
            -m 2G \
            -T results/rrna/{wildcards.sample}_tmp \
            -o {output.bam} \
            results/rrna/{wildcards.sample}.sam
        rm results/rrna/{wildcards.sample}.sam
        samtools flagstat {output.bam} > {output.flagstat}
        """


# =============================================================================
# STEP 6 — Kraken2 on unmapped reads
# =============================================================================
rule kraken2:
    input:
        r1 = "results/unmapped/{sample}_unmapped_R1.fq.gz",
        r2 = "results/unmapped/{sample}_unmapped_R2.fq.gz",
    output:
        report = "results/kraken2/{sample}.report",
        out    = "results/kraken2/{sample}.kraken",
    threads: 16
    resources:
        mem_mb  = 40000,
        runtime = 120,
    params:
        db = config["kraken2_db"],
    log:
        "logs/kraken2/{sample}.log"
    shell:
        """
        {ENV[kraken2]}
        kraken2 --paired \
            --db {params.db} \
            --threads {threads} \
            --report {output.report} \
            --gzip-compressed \
            {input.r1} {input.r2} \
            > {output.out} 2> {log}
        """


# =============================================================================
# STEP 7 — Trinity assembly per species (pooled unmapped replicates)
# =============================================================================
def get_species_unmapped(wc):
    samples = SAMPLES[wc.species]
    return {
        "r1": expand("results/unmapped/{sample}_unmapped_R1.fq.gz", sample=samples),
        "r2": expand("results/unmapped/{sample}_unmapped_R2.fq.gz", sample=samples),
    }


rule trinity_assembly:
    input:
        unpack(get_species_unmapped)
    output:
        fasta = "results/trinity/{species}_trinity.Trinity.fasta",
    threads: 16
    resources:
        mem_mb  = 120000,
        runtime = 1440,
    params:
        max_mem = "120G",
        outdir  = "results/trinity/{species}_trinity/",
    log:
        "logs/trinity/{species}.log"
    shell:
        """
        set +e  # disable exit on error for env sourcing
        {ENV[trinity]}
        set -e  # re-enable after env is loaded
        mkdir -p results/trinity/{wildcards.species}_trinity/
        cat {input.r1} > results/trinity/{wildcards.species}_trinity/all_R1.fq.gz
        cat {input.r2} > results/trinity/{wildcards.species}_trinity/all_R2.fq.gz
        Trinity --seqType fq \
            --left  results/trinity/{wildcards.species}_trinity/all_R2.fq.gz \
            --right results/trinity/{wildcards.species}_trinity/all_R1.fq.gz \
            --CPU {threads} \
            --max_memory {params.max_mem} \
            --normalize_reads \
            --output results/trinity/{wildcards.species}_trinity/ \
            2> {log}
        rm -f results/trinity/{wildcards.species}_trinity/all_R1.fq.gz \
              results/trinity/{wildcards.species}_trinity/all_R2.fq.gz
        """


rule trinity_stats:
    input:
        "results/trinity/{species}_trinity.Trinity.fasta"
    output:
        "results/trinity/{species}_trinity.assembly_stats.txt"
    log:
        "logs/trinity_stats/{species}.log"
    shell:
        """
        {ENV[trinity]}
        TrinityStats.pl {input} > {output} 2> {log}
        """


# =============================================================================
# STEP 8 — DIAMOND BLASTx (protein space — finds divergent viruses)
# =============================================================================
rule diamond_blastx:
    input:
        contigs = "results/trinity/{species}_trinity.Trinity.fasta",
        db      = config["diamond_db"],
    output:
        "results/diamond/{species}/blastx.tsv"
    threads: 16
    resources:
        mem_mb  = 32000,
        runtime = 480,
    params:
        evalue   = config.get("diamond_evalue",  "1e-5"),
        max_hits = config.get("diamond_max_hits", 10),
        extra    = config.get("diamond_extra",   "--more-sensitive"),
    log:
        "logs/diamond/{species}.log"
    shell:
        """
        {ENV[diamond]}
        diamond blastx \
            --query {input.contigs} \
            --db {input.db} \
            --out {output} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen \
                       qstart qend sstart send evalue bitscore stitle \
            --evalue {params.evalue} \
            --max-target-seqs {params.max_hits} \
            {params.extra} \
            --threads {threads} \
            2> {log}
        """


# =============================================================================
# STEP 9 — Taxonomy annotation
# =============================================================================
rule extract_accessions:
    input:
        diamond = "results/diamond/{species}/blastx.tsv",
    output:
        "results/taxonomy/{species}/subject_accessions.txt"
    shell:
        """
        awk '{{print $2}}' {input.diamond} \
            | sort -u \
            | sed 's/\\.[0-9]*$//' \
            > {output}
        """


rule map_accessions_to_taxid:
    input:
        accessions = "results/taxonomy/{species}/subject_accessions.txt",
        acc2taxid  = config["prot_acc2taxid"],
    output:
        "results/taxonomy/{species}/acc2taxid.tsv"
    threads: 4
    resources:
        mem_mb  = 16000,
        runtime = 60,
    log:
        "logs/taxonomy/{species}_acc2taxid.log"
    shell:
        """
        grep -F -f {input.accessions} {input.acc2taxid} \
            | awk '{{print $2"\\t"$3}}' \
            > {output} 2> {log}
        """


rule taxid_to_lineage:
    input:
        "results/taxonomy/{species}/acc2taxid.tsv"
    output:
        lineage = "results/taxonomy/{species}/lineages.tsv",
        summary = "results/taxonomy/{species}/taxonomy_summary.tsv",
    params:
        taxdump = config["taxdump_dir"],
    threads: 4
    resources:
        mem_mb  = 8000,
        runtime = 30,
    log:
        "logs/taxonomy/{species}_lineage.log"
    shell:
        """
        source activate snakemake_env
        awk '{{print $2}}' {input} | sort -u \
        | taxonkit lineage --data-dir {params.taxdump} \
        | taxonkit reformat --data-dir {params.taxdump} \
            -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" \
        > {output.lineage} 2> {log}

        join -t $'\\t' -1 2 -2 1 \
            <(sort -t$'\\t' -k2 {input}) \
            <(sort -t$'\\t' -k1 {output.lineage}) \
        | awk 'BEGIN{{OFS="\\t"; \
                print "Accession","TaxID","Kingdom","Phylum","Class","Order","Family","Genus","Species"}} \
               {{print $2,$1,$3,$4,$5,$6,$7,$8,$9}}' \
        > {output.summary} 2>> {log}
        """


# =============================================================================
# STEP 10 — Pipeline summary table
# =============================================================================
rule pipeline_summary:
    input:
        star_logs = expand("results/host_mapping/{sample}/Log.final.out", sample=ALL_SAMPLES),
        flagstats = expand("results/rrna/{sample}.flagstat",               sample=ALL_SAMPLES),
        kraken    = expand("results/kraken2/{sample}.report",              sample=ALL_SAMPLES),
        trinity   = expand("results/trinity/{species}_trinity.assembly_stats.txt", species=SPECIES),
        diamond   = expand("results/diamond/{species}/blastx.tsv",         species=SPECIES),
    output:
        "results/summary/pipeline_summary.tsv"
    run:
        import re, os

        header = ["sample", "species", "total_reads",
                  "pct_uniquely_mapped", "pct_unmapped",
                  "pct_rrna", "kraken_classified_pct"]
        rows = []

        for sample in ALL_SAMPLES:
            species = SAMPLE_TO_SPECIES[sample]
            row = {"sample": sample, "species": species}

            # STAR
            with open(f"results/host_mapping/{sample}/Log.final.out") as f:
                content = f.read()
            total  = re.search(r"Number of input reads \|\s+(\d+)",              content)
            mapped = re.search(r"Uniquely mapped reads % \|\s+([\d.]+)%",        content)
            unmap  = re.search(r"% of reads unmapped: too short \|\s+([\d.]+)%", content)
            row["total_reads"]         = total.group(1)  if total  else "NA"
            row["pct_uniquely_mapped"] = mapped.group(1) if mapped else "NA"
            row["pct_unmapped"]        = unmap.group(1)  if unmap  else "NA"

            # rRNA
            with open(f"results/rrna/{sample}.flagstat") as f:
                lines = f.readlines()
            total_r  = int(lines[0].split()[0]) if lines else 0
            mapped_r = int(lines[4].split()[0]) if len(lines) > 4 else 0
            row["pct_rrna"] = f"{100*mapped_r/total_r:.2f}" if total_r > 0 else "NA"

            # Kraken2
            with open(f"results/kraken2/{sample}.report") as f:
                first = f.readline()
            unclass = first.strip().split("\t")[0] if first else "NA"
            row["kraken_classified_pct"] = \
                f"{100 - float(unclass):.2f}" if unclass != "NA" else "NA"

            rows.append(row)

        os.makedirs("results/summary", exist_ok=True)
        with open(output[0], "w") as out:
            out.write("\t".join(header) + "\n")
            for row in rows:
                out.write("\t".join(str(row.get(h, "NA")) for h in header) + "\n")
# =============================================================================
# STEP 11 — Salmon quantification of viral contigs
# =============================================================================
rule salmon_index:
    input:
        "results/diamond/{species}/viral_contigs.fasta"
    output:
        directory("results/salmon/{species}/index/")
    threads: 8
    resources:
        mem_mb  = 16000,
        runtime = 30,
    log:
        "logs/salmon/{species}_index.log"
    shell:
        """
        {ENV[salmon_env]}
        salmon index \
            -t {input} \
            -i {output} \
            --threads {threads} \
            2> {log}
        """

rule salmon_quant:
    input:
        r1    = lambda wc: f"{config['raw_data_dir']}/{wc.sample}_R1_001.fastq.gz",
        r2    = lambda wc: f"{config['raw_data_dir']}/{wc.sample}_R2_001.fastq.gz",
        index = lambda wc: f"results/salmon/{SAMPLE_TO_SPECIES[wc.sample]}/index/",
    output:
        "results/salmon/{sample}/quant.sf"
    threads: 8
    resources:
        mem_mb  = 16000,
        runtime = 60,
    params:
        outdir = "results/salmon/{sample}/"
    log:
        "logs/salmon/{sample}.log"
    shell:
        """
        {ENV[salmon_env]}
        salmon quant \
            -i {input.index} \
            -l A \
            -1 {input.r1} \
            -2 {input.r2} \
            --threads {threads} \
            --validateMappings \
            -o {params.outdir} \
            2> {log}
        """

rule salmon_merge:
    input:
        lambda wc: expand("results/salmon/{sample}/quant.sf",
                          sample=SAMPLES[wc.species])
    output:
        "results/salmon/{species}/merged_tpm.tsv"
    run:
        import pandas as pd
        
        dfs = []
        for f in input:
            sample = f.split("/")[2]
            df = pd.read_csv(f, sep="\t", usecols=["Name","TPM","NumReads"])
            df = df.rename(columns={"TPM": f"{sample}_TPM",
                                    "NumReads": f"{sample}_reads"})
            dfs.append(df.set_index("Name"))
        
        merged = pd.concat(dfs, axis=1)
        merged.to_csv(output[0], sep="\t")
        print(f"Merged {len(input)} samples, {len(merged)} contigs")
        
rule merge_final_results:
    input:
        tpm = "results/salmon/{species}/merged_tpm.tsv",
        hits = "results/diamond/{species}/best_hits_filtered.tsv"
    output:
        "results/summary/{species}_viral_abundance_summary.tsv"
    run:
        import pandas as pd

        # Load the filtered DIAMOND hits
        hits_df = pd.read_csv(input.hits, sep="\t")

        # Load the Salmon TPM/Read counts
        # We use 'Name' as index because it corresponds to the Trinity 'qseqid'
        tpm_df = pd.read_csv(input.tpm, sep="\t").set_index("Name")

        # Merge the datasets on the contig ID
        # Assuming the ID column in best_hits_filtered.tsv is 'qseqid'
        final_df = hits_df.merge(tpm_df, left_on="qseqid", right_index=True, how="left")

        # Save the integrated table
        final_df.to_csv(output[0], sep="\t", index=False)
        print(f"Final summary generated for {wildcards.species}")