import os


def write_sample_map(log_path, individuals, gvcf_dir):
    """Write GenomicsDBImport sample-name-map: <sample>\\t<gvcf path>."""
    os.makedirs(log_path, exist_ok=True)
    out = os.path.join(log_path, 'sample_map.tsv')
    with open(out, 'w') as f:
        for ind in individuals:
            f.write(f'{ind}\t{gvcf_dir}/{ind}/{ind}.g.vcf.gz\n')
    return out
