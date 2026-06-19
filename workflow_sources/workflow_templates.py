from gwf import AnonymousTarget


def bwa2_index(ref, index_dir, basename, log_path, account):
    inputs = [ref]
    outputs = [f'{log_path}/{basename}_index.DONE']
    options = {'cores': 12, 'memory': '128g', 'walltime': '12:00:00', 'account': account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bwa2
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {index_dir} {log_path}
    cd {index_dir}
    bwa-mem2 index -p {basename} {ref}
    echo "DONE: $(date)"
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def gatk_dict(ref, log_path, account):
    dict_out = ref.rsplit('.', 1)[0] + '.dict'
    inputs = [ref]
    outputs = [f'{log_path}/gatk_dict.DONE']
    options = {'cores': 2, 'memory': '8g', 'walltime': '2:00:00', 'account': account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools117
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {log_path}
    [ -f {ref}.fai ] || samtools faidx {ref}
    conda activate gatk4
    if [ ! -f {dict_out} ]; then
        gatk CreateSequenceDictionary -R {ref} -O {dict_out}
    fi
    echo "DONE: $(date)"
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def bwa2_align(indname, r1, r2, index_prefix, out_dir, log_path, basename, account):
    inputs = [r1, r2, f'{log_path}/{basename}_index.DONE']
    outputs = [f'{log_path}/{indname}_align_s0.DONE']
    options = {'cores': 20, 'memory': '64g', 'walltime': '24:00:00', 'account': account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bwa2
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {out_dir} {log_path}
    cd {out_dir}
    bwa-mem2 mem -t 20 {index_prefix} {r1} {r2} | samtools view -Sb -@ 4 - > {indname}_s0.bam
    echo "DONE: $(date)"
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def picard_addRG(indname, out_dir, log_path, account):
    inputs = [f'{log_path}/{indname}_align_s0.DONE']
    outputs = [f'{log_path}/{indname}_addRG_s1.DONE']
    options = {'cores': 4, 'memory': '16g', 'walltime': '12:00:00', 'account': account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate gatk4
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    cd {out_dir}
    gatk AddOrReplaceReadGroups -I {indname}_s0.bam -O {indname}_s1.bam \\
        -RGID {indname} -RGPU unknown -RGSM {indname} -RGPL illumina -RGLB lib0
    rm {indname}_s0.bam
    echo "DONE: $(date)"
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def samtools_fixmate(indname, out_dir, log_path, account):
    inputs = [f'{log_path}/{indname}_addRG_s1.DONE']
    outputs = [f'{log_path}/{indname}_fixmate_s2.DONE']
    options = {'cores': 16, 'memory': '32g', 'walltime': '12:00:00', 'account': account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools117
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    cd {out_dir}
    samtools fixmate -rm -@ 16 {indname}_s1.bam {indname}_s2.bam
    rm {indname}_s1.bam
    echo "DONE: $(date)"
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def samtools_sort(indname, out_dir, log_path, account):
    inputs = [f'{log_path}/{indname}_fixmate_s2.DONE']
    outputs = [f'{log_path}/{indname}_sort_s3.DONE']
    options = {'cores': 16, 'memory': '32g', 'walltime': '12:00:00', 'account': account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools117
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    cd {out_dir}
    samtools sort -@ 16 -o {indname}_s3.bam {indname}_s2.bam
    rm {indname}_s2.bam
    echo "DONE: $(date)"
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def samtools_markdup(indname, out_dir, log_path, account):
    inputs = [f'{log_path}/{indname}_sort_s3.DONE']
    outputs = [f'{log_path}/{indname}_markdup_s4.DONE']
    options = {'cores': 16, 'memory': '32g', 'walltime': '12:00:00', 'account': account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools117
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    cd {out_dir}
    samtools markdup -r -f {indname}.markdup.stat -s -@ 16 {indname}_s3.bam {indname}_s4.bam
    rm {indname}_s3.bam
    echo "DONE: $(date)"
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def samtools_mq_filter(indname, out_dir, log_path, account):
    inputs = [f'{log_path}/{indname}_markdup_s4.DONE']
    outputs = [f'{log_path}/{indname}_MQ_s5.DONE']
    options = {'cores': 16, 'memory': '32g', 'walltime': '12:00:00', 'account': account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools117
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    cd {out_dir}
    samtools view -@ 16 -bq 60 -f 0x2 -F 0x4 {indname}_s4.bam > {indname}_final.bam
    rm {indname}_s4.bam
    echo "DONE: $(date)"
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def samtools_index(indname, out_dir, log_path, account):
    inputs = [f'{log_path}/{indname}_MQ_s5.DONE']
    outputs = [f'{log_path}/{indname}_bai_s6.DONE']
    options = {'cores': 8, 'memory': '8g', 'walltime': '4:00:00', 'account': account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools117
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    cd {out_dir}
    samtools index -@ 8 {indname}_final.bam
    echo "DONE: $(date)"
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def gatk_haplotype_call(indname, tag, interval, ref, bam_path, gvcf_dir, log_path, account):
    inputs = [
        f'{log_path}/{indname}_bai_s6.DONE',
        f'{log_path}/gatk_dict.DONE',
    ]
    outputs = [f'{log_path}/{indname}_{tag}_gvcf.DONE']
    options = {'cores': 4, 'memory': '24g', 'walltime': '36:00:00', 'account': account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate gatk4
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {gvcf_dir}
    cd {gvcf_dir}
    gatk --java-options "-Xmx20g" HaplotypeCaller \\
        -R {ref} \\
        -I {bam_path} \\
        -O {indname}_{tag}.g.vcf.gz \\
        -L {interval} \\
        -ERC GVCF \\
        --tmp-dir /scratch/$SLURM_JOBID/ \\
        --native-pair-hmm-threads 4
    echo "DONE: $(date)"
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def merge_gvcf_by_chrom(indname, tags, gvcf_dir, log_path, account):
    inputs = [f'{log_path}/{indname}_{c}_gvcf.DONE' for c in tags]
    outputs = [f'{log_path}/{indname}_gvcf_merged.DONE']
    options = {'cores': 6, 'memory': '24g', 'walltime': '12:00:00', 'account': account}
    concat = ' '.join(f'{gvcf_dir}/{indname}_{c}.g.vcf.gz' for c in tags)
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bcftools
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    cd {gvcf_dir}
    bcftools concat --threads 6 -O z -o {indname}.g.vcf.gz {concat}
    conda activate gatk4
    gatk --java-options "-Xmx20g" IndexFeatureFile --input {indname}.g.vcf.gz
    echo "DONE: $(date)"
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def gatk_consolidate(tag, interval, individuals, ref, vcf_dir, sample_map, log_path, batch_size, account):
    inputs = [f'{log_path}/{ind}_gvcf_merged.DONE' for ind in individuals]
    outputs = [f'{log_path}/joint_{tag}.DONE']
    options = {'cores': 6, 'memory': '200g', 'walltime': '24:00:00', 'account': account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate gatk4
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {vcf_dir}
    gatk --java-options "-Xmx60g -Xms60g" GenomicsDBImport \\
        --genomicsdb-workspace-path /scratch/$SLURM_JOBID/db_{tag} \\
        --batch-size {batch_size} \\
        --sample-name-map {sample_map} \\
        --tmp-dir /scratch/$SLURM_JOBID/ \\
        --reader-threads 6 \\
        -L {interval} \\
        --genomicsdb-shared-posixfs-optimizations true \\
        --genomicsdb-vcf-buffer-size 4194304
    gatk --java-options "-Xmx60g -Xms60g" GenotypeGVCFs \\
        -R {ref} \\
        -V gendb:///scratch/$SLURM_JOBID/db_{tag} \\
        --tmp-dir /scratch/$SLURM_JOBID/ \\
        -O /scratch/$SLURM_JOBID/{tag}.vcf.gz \\
        -L {interval} \\
        -all-sites true
    gatk --java-options "-Xmx60g -Xms60g" VariantFiltration \\
        -V /scratch/$SLURM_JOBID/{tag}.vcf.gz \\
        -O {vcf_dir}/{tag}_GATK_filtered.vcf.gz \\
        --filter-name "gatk_germline" \\
        --tmp-dir /scratch/$SLURM_JOBID/ \\
        --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \\
        --verbosity ERROR
    echo "DONE: $(date)"
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def merge_vcfs(tags, vcf_dir, joint_name, log_path, account):
    inputs = [f'{log_path}/joint_{c}.DONE' for c in tags]
    outputs = [f'{log_path}/{joint_name}_concat.DONE']
    options = {'cores': 6, 'memory': '24g', 'walltime': '12:00:00', 'account': account}
    concat = ' '.join(f'{vcf_dir}/{c}_GATK_filtered.vcf.gz' for c in tags)
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bcftools
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    cd {vcf_dir}
    bcftools concat --threads 6 -O z -o {joint_name}.vcf.gz {concat}
    conda activate gatk4
    gatk --java-options "-Xmx20g" IndexFeatureFile --input {joint_name}.vcf.gz
    echo "DONE: $(date)"
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
