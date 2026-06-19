import os
import yaml

from workflow_templates import (
    bwa2_index,
    gatk_dict,
    bwa2_align,
    picard_addRG,
    samtools_fixmate,
    samtools_sort,
    samtools_markdup,
    samtools_mq_filter,
    samtools_index,
    gatk_haplotype_call,
    merge_gvcf_by_chrom,
    gatk_consolidate,
    merge_vcfs,
)
from sample_map_writer import write_sample_map


def _read_contigs(fai_path):
    """Return an ordered list of (contig_name, length) from a .fai index."""
    contigs = []
    with open(fai_path) as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.split('\t')
            contigs.append((fields[0], int(fields[1])))
    return contigs


def _make_units(fai_path, window_size):
    """Build the list of calling units (genomic intervals to parallelize over).

    Each unit is a (tag, interval, name_safe) tuple where:
      - tag       : token used in output filenames / .DONE sentinels.
      - interval  : value passed to GATK -L.
      - name_safe : token used in gwf target names ('.' is not target-name safe).

    window_size <= 0 -> one unit per whole contig (the default, contig-level
    parallelization). window_size > 0 -> each contig is split into consecutive
    non-overlapping windows of at most window_size bp (1-based, inclusive),
    yielding many small units for finer ("massive") parallelization.
    """
    units = []
    for contig, length in _read_contigs(fai_path):
        if window_size and window_size > 0:
            start = 1
            while start <= length:
                end = min(start + window_size - 1, length)
                tag = f'{contig}_{start}_{end}'
                units.append((tag, f'{contig}:{start}-{end}', tag.replace('.', '_')))
                start = end + 1
        else:
            units.append((contig, contig, contig.replace('.', '_')))
    return units


def _read_sample_sheet(path):
    """Parse a tab-separated sample sheet.

    Expected columns (tab-delimited): sample_id, r1, r2

    The first non-comment line may be a header (its first field equal to
    'sample_id', case-insensitive); it is skipped if present. Blank lines and
    lines starting with '#' are ignored.

    Returns an ordered list of (sample_id, r1, r2) tuples.
    """
    samples = []
    seen = set()
    with open(path) as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.rstrip('\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            fields = [c.strip() for c in line.split('\t')]
            if fields[0].lower() == 'sample_id':
                continue
            if len(fields) < 3 or not fields[0] or not fields[1] or not fields[2]:
                raise ValueError(
                    f'{path}:{lineno}: expected 3 tab-separated columns '
                    f'(sample_id, r1, r2), got: {line!r}'
                )
            sample_id, r1, r2 = fields[0], fields[1], fields[2]
            if sample_id in seen:
                raise ValueError(f'{path}:{lineno}: duplicate sample_id {sample_id!r}')
            seen.add(sample_id)
            samples.append((sample_id, r1, r2))
    if not samples:
        raise ValueError(f'{path}: no samples found')
    return samples


def gatk_calling_workflow(config_file, gwf):
    with open(config_file) as f:
        cfg = yaml.safe_load(f)

    account = cfg['account']
    log_path = cfg['log_directory_path']
    out_root = cfg['output_directory_path']
    ref = cfg['reference_genome']
    sample_sheet = cfg['sample_sheet']
    basename = cfg['index_basename']
    joint_name = cfg['joint_name']
    batch_size = cfg['batch_size']
    window_size = cfg.get('window_size', 0)

    samples = _read_sample_sheet(sample_sheet)
    individuals = [s[0] for s in samples]

    units = _make_units(cfg['reference_fai'], window_size)
    tags = [u[0] for u in units]

    index_dir = os.path.join(out_root, 'index')
    index_prefix = os.path.join(index_dir, basename)
    align_root = os.path.join(out_root, 'align')
    gvcf_root = os.path.join(out_root, 'gvcf')
    joint_dir = os.path.join(out_root, 'joint')

    os.makedirs(log_path, exist_ok=True)

    gwf.target_from_template(
        name='index_ref',
        template=bwa2_index(ref, index_dir, basename, log_path, account),
    )
    gwf.target_from_template(
        name='gatk_dict',
        template=gatk_dict(ref, log_path, account),
    )

    for ind, r1, r2 in samples:
        ind_align_dir = os.path.join(align_root, ind)
        ind_gvcf_dir = os.path.join(gvcf_root, ind)
        final_bam = os.path.join(ind_align_dir, f'{ind}_final.bam')

        gwf.target_from_template(
            name=f'{ind}_align',
            template=bwa2_align(ind, r1, r2, index_prefix, ind_align_dir, log_path, basename, account),
        )
        gwf.target_from_template(
            name=f'{ind}_addRG',
            template=picard_addRG(ind, ind_align_dir, log_path, account),
        )
        gwf.target_from_template(
            name=f'{ind}_fixmate',
            template=samtools_fixmate(ind, ind_align_dir, log_path, account),
        )
        gwf.target_from_template(
            name=f'{ind}_sort',
            template=samtools_sort(ind, ind_align_dir, log_path, account),
        )
        gwf.target_from_template(
            name=f'{ind}_markdup',
            template=samtools_markdup(ind, ind_align_dir, log_path, account),
        )
        gwf.target_from_template(
            name=f'{ind}_mq',
            template=samtools_mq_filter(ind, ind_align_dir, log_path, account),
        )
        gwf.target_from_template(
            name=f'{ind}_bai',
            template=samtools_index(ind, ind_align_dir, log_path, account),
        )
        for tag, interval, name_safe in units:
            gwf.target_from_template(
                name=f'{ind}_{name_safe}_gvcf',
                template=gatk_haplotype_call(
                    ind, tag, interval, ref, final_bam, ind_gvcf_dir, log_path, account
                ),
            )
        gwf.target_from_template(
            name=f'{ind}_gvcf_merge',
            template=merge_gvcf_by_chrom(ind, tags, ind_gvcf_dir, log_path, account),
        )

    sample_map = write_sample_map(log_path, individuals, gvcf_root)

    for tag, interval, name_safe in units:
        gwf.target_from_template(
            name=f'joint_{name_safe}',
            template=gatk_consolidate(
                tag, interval, individuals, ref, joint_dir, sample_map, log_path, batch_size, account
            ),
        )

    gwf.target_from_template(
        name=f'{joint_name}_concat',
        template=merge_vcfs(tags, joint_dir, joint_name, log_path, account),
    )

    return gwf
