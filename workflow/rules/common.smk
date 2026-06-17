# common.smk — shared helpers used by both modes.

wildcard_constraints:
    sample=r"[A-Za-z0-9_.\-]+",


def sample_kind(sample):
    return SAMPLES[sample]["kind"]


def sample_path(sample):
    return SAMPLES[sample]["path"]


def sample_path2(sample):
    return SAMPLES[sample].get("path2", "")


def sample_is_protein_mode(sample):
    return sample_kind(sample) in ("protein", "nucleotide")


def sample_is_read_mode(sample):
    return sample_kind(sample) == "fastq"


# Used by report.smk to dispatch the right calls-source per sample.
def calls_tsv_for(sample):
    return str(RESULTS / sample / "calls" / "scycle_calls.tsv")
