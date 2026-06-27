import boto3
import os
import subprocess
import json
import re
import threading
import shlex
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import pandas as pd

# NOTE: emailing has been removed. The web app addresses runs purely by their
# result code (the flat S3 slug), and the results viewer polls the bucket, so
# there is no send_email import or call anywhere in this pipeline anymore.
from gen_presign_url import generate_presigned_url, shorten_url  # noqa: F401  (kept for optional CLI use)
import mimetypes

# Initialize S3 client
s3_client = boto3.client('s3')

# Global list to track downloaded folders
downloaded_folders = []

MAX_PARALLEL_JOBS = 8
BRESEQ_THREADS_PER_JOB = 6
log_lock = threading.Lock()

log_file_path = '/home/ark/MAB/breseq/processed_folders.log'
failed_log_file_path = '/home/ark/MAB/breseq/failed_folders.log'

# Two-bucket layout:
#   INPUT_BUCKET   - the web app uploads run inputs here (reads, reference,
#                    form-data.txt) under a flat <slug>/ prefix. The tower
#                    LISTS and DOWNLOADS from here.
#   RESULTS_BUCKET - the tower WRITES all results here under the SAME <slug>/
#                    prefix. The web app's results viewer reads from here
#                    (it needs public ListBucket + GetObject + CORS).
# Both can be overridden via environment variables.
INPUT_BUCKET = os.environ.get('BRESEQ_INPUT_BUCKET', 'midauthorbio-breseq-input')
RESULTS_BUCKET = os.environ.get('BRESEQ_RESULTS_BUCKET', 'midauthorbio-breseq-results')

speciesDict = {"Pseudomonas fluorescens SBW25": "GCA_931907645.1",
                "Escherichia coli K-12 MG1655": "GCA_000005845.2",
               "Bacillus subtilis 168": "GCA_000009045.1",
               "Pseudomonas aeruginosa PAO1": "GCA_000006765.1"}

app = "breseq"

preset_dir = "/home/ark/MAB/breseq/fastq_presets"

# breseq flags we accept from the web app's `BreseqFlags:` line. The frontend
# only ever emits flags from a fixed catalog, but because this string is spliced
# into the command we re-validate every token against this allow-list and drop
# anything unexpected. Flags that take a value are followed by a single token.
ALLOWED_BRESEQ_FLAGS = {
    # value-taking (numeric) flags
    "--minimum-mapping-quality": 1,
    "--base-quality-cutoff": 1,
    "--require-match-length": 1,
    "--require-match-fraction": 1,
    "--maximum-read-mismatches": 1,
    "--deletion-coverage-propagation-cutoff": 1,
    "--deletion-coverage-seed-cutoff": 1,
    "--junction-indel-split-length": 1,
    "--junction-alignment-pair-limit": 1,
    "--junction-minimum-candidates": 1,
    "--junction-maximum-candidates": 1,
    "--junction-candidate-length-factor": 1,
    "--junction-minimum-candidate-pos-hash-score": 1,
    "--junction-score-cutoff": 1,
    "--junction-minimum-pos-hash-score": 1,
    "--junction-minimum-side-match": 1,
    "--junction-minimum-pr-no-read-start-per-position": 1,
    "--consensus-score-cutoff": 1,
    "--consensus-frequency-cutoff": 1,
    "--consensus-minimum-variant-coverage": 1,
    "--consensus-minimum-total-coverage": 1,
    "--consensus-minimum-variant-coverage-each-strand": 1,
    "--consensus-minimum-total-coverage-each-strand": 1,
    "--consensus-reject-indel-homopolymer-length": 1,
    "--consensus-reject-surrounding-homopolymer-length": 1,
    "--polymorphism-score-cutoff": 1,
    "--polymorphism-frequency-cutoff": 1,
    "--polymorphism-minimum-variant-coverage": 1,
    "--polymorphism-minimum-total-coverage": 1,
    "--polymorphism-minimum-variant-coverage-each-strand": 1,
    "--polymorphism-minimum-total-coverage-each-strand": 1,
    "--polymorphism-bias-cutoff": 1,
    "--polymorphism-reject-indel-homopolymer-length": 1,
    "--polymorphism-reject-surrounding-homopolymer-length": 1,
    # boolean flags (no value)
    "--quality-score-trim": 0,
    "--no-junction-prediction": 0,
    "--junction-allow-suboptimal-matches": 0,
    "--polymorphism-no-indels": 0,
}

_NUMERIC_RX = re.compile(r"^-?\d+(\.\d+)?$")


def validate_fastq_file(fastq_path):
    """
    Basic FASTQ integrity check.

    Returns:
        (True, None) if valid
        (False, error_message) if invalid
    """

    if not fastq_path:
        return False, "FASTQ path is missing"

    if not os.path.exists(fastq_path):
        return False, f"FASTQ file does not exist: {fastq_path}"

    if os.path.getsize(fastq_path) == 0:
        return False, f"FASTQ file is empty: {fastq_path}"

    # For plain FASTQ only. gzipped FASTQs are not line-count checked here.
    if fastq_path.endswith(".gz"):
        return True, None

    with open(fastq_path, "r", errors="replace") as f:
        line_count = sum(1 for _ in f)

    if line_count % 4 != 0:
        return False, (
            f"FASTQ appears incomplete/corrupt: {fastq_path}. "
            f"Line count is {line_count}, which is not divisible by 4."
        )

    return True, None


def _kv(line):
    """Split a 'Key: value' form-data line into (key, value).

    The web app writes 'Key: value' (colon-delimited). We split on the FIRST
    colon so values that contain spaces or extra colons survive intact.
    Returns (None, None) for lines without a colon.
    """
    if ":" not in line:
        return None, None
    key, value = line.split(":", 1)
    return key.strip(), value.strip()


def _sanitize_flags(raw_flags):
    """Validate the web app's BreseqFlags string against the allow-list.

    Returns a clean list of tokens (e.g. ['--minimum-mapping-quality', '20',
    '--no-junction-prediction']). Unknown flags or malformed values are dropped
    with a warning so a bad/hand-edited form-data.txt can't inject arbitrary
    arguments.
    """
    tokens = shlex.split(raw_flags) if raw_flags else []
    clean = []
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok not in ALLOWED_BRESEQ_FLAGS:
            print(f"[flags] Dropping unrecognized token: {tok!r}")
            i += 1
            continue
        takes_value = ALLOWED_BRESEQ_FLAGS[tok]
        if takes_value:
            if i + 1 >= len(tokens):
                print(f"[flags] {tok} expects a value but none given; dropping")
                break
            val = tokens[i + 1]
            if not _NUMERIC_RX.match(val):
                print(f"[flags] {tok} got non-numeric value {val!r}; dropping")
                i += 2
                continue
            clean.extend([tok, val])
            i += 2
        else:
            clean.append(tok)
            i += 1
    return clean


def extract_form_data(folder_path):
    """Parse the web app's form-data.txt.

    New contract (no email, flat slug, explicit samples, advanced flags):

        Accession: GCA_000005845.2            # or N/A
        ReferenceFile: REL606.gbk             # or N/A
        Sample: clone1 | R1=clone1_R1.fastq.gz | R2=clone1_R2.fastq.gz
        Sample: clone2 | R1=clone2_R1.fastq.gz | R2=N/A
        ForwardReads: clone1_R1.fastq.gz,clone2_R1.fastq.gz     # flattened
        ReverseReads: clone1_R2.fastq.gz,N/A                    # flattened or N/A
        Polymorphic: clonal                   # clonal | population
        Name: REL606_clone1                   # -> breseq -n <name>
        ReferenceFormat: genbank              # informational
        BreseqFlags: --minimum-mapping-quality 20 --no-junction-prediction

    Returns:
        referenceFile : str   absolute path to the reference (or "None")
        poly          : str   "clonal" | "population"
        read_paths    : list  absolute paths to every read file, in order
        name          : str   run name for breseq -n (may be "")
        extra_flags   : list  validated advanced breseq flags
        app           : str   "breseq" (kept for compatibility)
        sample_specs  : list  one dict per sample, each:
                                {"label": <display label>,
                                 "key":   <fs/url-safe label>,
                                 "reads": [<abs r1>, <abs r2?>]}
                              Each sample is run as its OWN breseq job so the
                              results viewer can toggle between samples.
    """
    app = "breseq"
    form_file = os.path.join(folder_path, "form-data.txt")

    reference = accession = species = None
    poly = "clonal"
    name = ""
    raw_flags = ""
    samples = []          # list of (label, r1, r2_or_None)
    fwd_list = []         # flattened fallback
    rev_list = []

    if os.path.exists(form_file):
        with open(form_file, "r") as f:
            for raw_line in f:
                line = raw_line.rstrip("\n")
                key, value = _kv(line)
                if key is None:
                    continue

                if key == "Accession":
                    accession = value
                elif key == "ReferenceFile":
                    reference = value
                elif key == "SpeciesPreset":
                    species = value
                elif key == "Polymorphic":
                    poly = value
                elif key == "Name":
                    name = value
                elif key == "BreseqFlags":
                    raw_flags = value
                elif key == "Sample":
                    # "<label> | R1=<f> | R2=<f|N/A>"
                    label, r1, r2 = _parse_sample_value(value)
                    if r1:
                        samples.append((label, r1, r2))
                elif key == "ForwardReads":
                    fwd_list = [x.strip() for x in value.split(",")
                                if x.strip() and x.strip() != "N/A"]
                elif key == "ReverseReads":
                    rev_list = [x.strip() for x in value.split(",")
                                if x.strip() and x.strip() != "N/A"]
                elif key == "evolving":
                    app = "evolvingstem"

    # ---- Build the ordered list of read files (R1, R2, R1, R2, ...) --------
    read_names = []
    if samples:
        for _label, r1, r2 in samples:
            read_names.append(r1)
            if r2 and r2 != "N/A":
                read_names.append(r2)
    else:
        # Fall back to the flattened lists (back-compat / single sample).
        # Interleave so each forward is followed by its matching reverse.
        for idx, r1 in enumerate(fwd_list):
            read_names.append(r1)
            if idx < len(rev_list):
                read_names.append(rev_list[idx])
        # any extra reverse files beyond the forwards
        for r2 in rev_list[len(fwd_list):]:
            read_names.append(r2)

    read_paths = [os.path.join(folder_path, n) for n in read_names]

    # ---- Build per-sample specs (each runs as its own breseq job) ----------
    # Prefer the explicit Sample: lines. If none were given, treat the whole
    # submission as a single sample (back-compat with older single-sample runs).
    sample_specs = []
    if samples:
        used_keys = set()
        for idx, (label, r1, r2) in enumerate(samples):
            disp = (label or f"sample{idx + 1}").strip()
            key = _safe_sample_key(disp, idx, used_keys)
            reads = [os.path.join(folder_path, r1)]
            if r2 and r2 != "N/A":
                reads.append(os.path.join(folder_path, r2))
            sample_specs.append({"label": disp, "key": key, "reads": reads})
    elif read_paths:
        disp = (name or "sample1").strip() or "sample1"
        sample_specs.append({
            "label": disp,
            "key": _safe_sample_key(disp, 0, set()),
            "reads": read_paths,
        })

    # ---- Resolve the reference ---------------------------------------------
    if reference and reference != "N/A":
        referenceFile = os.path.join(folder_path, reference)
    elif accession and accession != "N/A":
        os.system(f"/home/ark/MAB/bin/breseq-local/bit2local.sh -a {accession} -o {folder_path}")
        reference = accession + ".gb"
        referenceFile = os.path.join(folder_path, reference)
    elif species and species != "N/A" and species in speciesDict:
        speciesAcc = speciesDict[species]
        referenceFile = f"/home/ark/MAB/breseq/references/{speciesAcc}.gb"
    else:
        referenceFile = "None"

    extra_flags = _sanitize_flags(raw_flags)

    return referenceFile, poly, read_paths, name, extra_flags, app, sample_specs


def _safe_sample_key(label, idx, used_keys):
    """Turn a sample label into a filesystem/URL-safe, unique sub-prefix key.

    e.g. "A 2019/1" -> "A_2019_1". Falls back to sample<idx+1> and de-dupes.
    """
    key = re.sub(r"[^A-Za-z0-9._-]+", "_", (label or "").strip()).strip("_")
    if not key:
        key = f"sample{idx + 1}"
    base = key
    n = 2
    while key in used_keys:
        key = f"{base}_{n}"
        n += 1
    used_keys.add(key)
    return key


def _parse_sample_value(value):
    """Parse 'label | R1=<f> | R2=<f|N/A>' -> (label, r1, r2_or_None)."""
    label, r1, r2 = "", None, None
    parts = [p.strip() for p in value.split("|")]
    for i, part in enumerate(parts):
        if part.upper().startswith("R1="):
            r1 = part.split("=", 1)[1].strip()
        elif part.upper().startswith("R2="):
            v = part.split("=", 1)[1].strip()
            r2 = None if v == "N/A" else v
        elif i == 0:
            label = part
    return label, r1, r2


def generate_mutation_json(output_dir):
    """
    Generate mutation_predictions.json from breseq output.gd
    """

    output_dir = Path(output_dir)
    gd_file = output_dir / "output" / "output.gd"
    json_file = output_dir / "mutation_predictions.json"

    if not gd_file.exists():
        print(f"[mutation-json] ERROR: {gd_file} not found")
        return False

    mutations = []

    with open(gd_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            record_type = parts[0]

            # Only real mutation records
            if record_type not in {"SNP", "DEL", "INS", "SUB", "AMP", "MOB", "INV", "CON"}:
                continue

            try:
                seq_id = parts[2]
                position = int(parts[3])
            except (IndexError, ValueError):
                continue

            # Parse key=value fields
            info = {}
            for field in parts[4:]:
                if "=" in field:
                    k, v = field.split("=", 1)
                    info[k] = v

            mutations.append({
                "seq_id": seq_id,
                "position": position,
                "mutation": info.get("mutation", record_type),
                "annotation": info.get("annotation", ""),
                "gene": info.get("gene", ""),
                "description": info.get("product", "")
            })

    with open(json_file, "w") as f:
        json.dump(mutations, f, indent=2)

    print(f"[mutation-json] Wrote {len(mutations)} mutations → {json_file}")
    return True


def load_seen_folders(log_path):
    if os.path.exists(log_path):
        with open(log_path, 'r') as f:
            return set(line.strip() for line in f)
    return set()

def append_seen_folder(log_path, folder):
    with log_lock:
        with open(log_path, 'a') as f:
            f.write(folder + '\n')

def list_folders_in_bucket(bucket_name):
    """Find every run folder by locating <slug>/form-data.txt at the bucket root.

    The web app now uploads to a FLAT slug (no <user>/<subfolder> nesting), so a
    run key looks like 'AbCd1234Ef/form-data.txt'. We return 'AbCd1234Ef/'.
    """
    paginator = s3_client.get_paginator('list_objects_v2')
    response_iterator = paginator.paginate(Bucket=bucket_name)

    folders = set()
    for page in response_iterator:
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.endswith("form-data.txt"):
                parts = key.split('/')
                # flat layout: ["<slug>", "form-data.txt"]
                if len(parts) == 2 and parts[0]:
                    folders.add(f"{parts[0]}/")
    return sorted(folders)


def download_s3_folder(bucket_name, s3_folder, local_folder):
    os.makedirs(local_folder, exist_ok=True)

    paginator = s3_client.get_paginator("list_objects_v2")

    for page in paginator.paginate(Bucket=bucket_name, Prefix=s3_folder):
        for obj in page.get("Contents", []):
            key = obj["Key"]

            # SKIP DIRECTORY MARKERS AND PREFIX ITSELF
            if key == s3_folder or key == f"{s3_folder}/" or key.endswith("/"):
                continue

            # Do NOT re-download previously generated results (output/, tarball,
            # mutation json, coverage, etc.). Only the original inputs live at
            # the top level of the slug; results sit under output/ or are named
            # artifacts. We keep this simple: skip the output/ subtree.
            relative_path = os.path.relpath(key, s3_folder)
            if relative_path in (".", ""):
                continue
            if relative_path.startswith("output/"):
                continue

            local_file_path = os.path.join(local_folder, relative_path)
            os.makedirs(os.path.dirname(local_file_path), exist_ok=True)

            print(f"Downloading s3://{bucket_name}/{key} → {local_file_path}")

            s3_client.download_file(bucket_name, key, local_file_path)


def find_fastq_files(folder_path):
    fastq_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path)
                if f.endswith('.fastq') or f.endswith('.fastq.gz')]
    return fastq_files


def run_breseq_command(folder_path, read_paths, output_dir, poly, gbk_file,
                       name="", extra_flags=None, threads=BRESEQ_THREADS_PER_JOB):
    """
    Run breseq with one or more read files (single- or paired-end, multi-sample).

    read_paths : ordered list of absolute FASTQ paths (R1, R2, R1, R2, ...)
    name       : breseq -n run name (optional)
    extra_flags: validated advanced flag tokens (optional)
    """
    extra_flags = extra_flags or []
    os.makedirs(output_dir, exist_ok=True)

    fastq_files = [p for p in read_paths if p and p != "N/A"]
    if not fastq_files:
        raise FileNotFoundError("No read files provided to breseq")

    # ---- sanity checks (fail early, not later) ----
    for f in fastq_files + [gbk_file]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Missing required file: {f}")

    # ---- build breseq command SAFELY (LIST, not string) ----
    cmd = [
        "breseq",
        "-j", str(threads),
        "-l", "70",
        "-o", output_dir,
        "-r", gbk_file,
    ]

    # polymorphism mode (-p) when the sample is not clonal
    if poly != "clonal":
        cmd.append("--polymorphism-prediction")

    # run name (-n)
    if name:
        cmd.extend(["-n", name])

    # advanced, pre-validated flags from the web app
    cmd.extend(extra_flags)

    # add FASTQs (all of them, positional)
    cmd.extend(fastq_files)

    print("[breseq] Running:")
    print(" ".join(cmd))

    # ---- run inside conda env WITHOUT bash -c ----
    full_command = [
        "/home/ark/miniconda3/bin/conda",
        "run",
        "-n",
        "breseq_env",
        *cmd
    ]

    result = subprocess.run(full_command, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"[breseq] ERROR in {folder_path}")
        print(result.stderr)
        return False

    # ---- verify output actually exists (prevents cascade failures) ----
    gd_file = os.path.join(output_dir, "output", "output.gd")
    if not os.path.exists(gd_file):
        print(f"[breseq] FAILED — missing output.gd in {output_dir}")
        return False

    mode = "paired-end" if len(fastq_files) >= 2 else "single-end"
    print(f"[breseq] Success ({mode}, {len(fastq_files)} file(s)) → {output_dir}")
    return True



def run_samtools_command(output_dir):
    bam_file = os.path.join(output_dir, "data", "reference.bam")
    coverage_file = os.path.join(output_dir, "data", "coverage.txt")

    if not os.path.exists(bam_file):
        print(f"Error: reference.bam not found in {bam_file}")
        return None

    # Run samtools depth command
    command = f"samtools depth -a {bam_file} > {coverage_file}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running samtools: {result.stderr}")
        return None
    else:
        print(f"Samtools depth ran successfully.")
    return coverage_file


def calculate_coverage_averages(coverage_file, output_dir):
    averages_file = os.path.join(output_dir, "averages.csv")

    if not os.path.exists(coverage_file):
        print(f"Error: coverage.txt not found in {coverage_file}")
        return

    df = pd.read_csv(coverage_file, sep="\t", header=None, names=["ID", "Index", "Value"])
    averages = []
    chunk_size = 1000

    for i in range(0, len(df), chunk_size):
        chunk = df.iloc[i:i + chunk_size]
        avg = round(chunk["Value"].mean(), 2)
        averages.append((i + 1, i + len(chunk), avg))

    averages_df = pd.DataFrame(averages, columns=["Start Row", "End Row", "Average Value"])
    averages_df.to_csv(averages_file, index=False)
    print(f"Coverage averages saved to {averages_file}")


def upload_file_to_s3(bucket_name, s3_folder, local_file):
    s3_key = os.path.join(s3_folder, os.path.basename(local_file))
    s3_client.upload_file(local_file, bucket_name, s3_key)
    print(f"Uploaded {local_file} to s3://{bucket_name}/{s3_key}")

def upload_directory_with_mime(local_dir, bucket, prefix):

    s3 = boto3.client('s3')
    import mimetypes
    import os

    if not os.path.isdir(local_dir):
        raise RuntimeError(f"HTML output directory does not exist: {local_dir}")

    for root, _, files in os.walk(local_dir):
        for filename in files:
            file_path = os.path.join(root, filename)

            # RELATIVE path inside output/
            relative_path = os.path.relpath(file_path, local_dir)

            # FINAL S3 KEY (NO LOCAL PATHS)
            s3_key = f"{prefix}/{relative_path}".replace("\\", "/")

            content_type, _ = mimetypes.guess_type(file_path)
            if content_type is None:
                content_type = "binary/octet-stream"

            print(f"Uploading {file_path} → s3://{bucket}/{s3_key}")

            s3.upload_file(
                Filename=file_path,
                Bucket=bucket,
                Key=s3_key,
                ExtraArgs={
                    "ContentType": content_type,
                    "ContentDisposition": "inline"
                }
            )

def s3_key_exists(bucket_name, key):
    try:
        s3_client.head_object(Bucket=bucket_name, Key=key)
        return True
    except s3_client.exceptions.ClientError as e:
        error_code = e.response.get("Error", {}).get("Code", "")
        if error_code in ("404", "NoSuchKey", "NotFound"):
            return False
        raise


def _process_one_sample(sample, s3_folder, results_bucket, referenceFile, poly,
                        name, extra_flags, local_base_dir, base_output_dir):
    """Run breseq for ONE sample and upload its results to a per-sample sub-prefix.

    Results land at:
        <slug>/<sample_key>/output/index.html   (+ summary.html, marginal.html, ...)
        <slug>/<sample_key>/mutation_predictions.json
        <slug>/<sample_key>/<sample_key>.tar.gz
        <slug>/<sample_key>/coverage.txt, averages.csv, reference.bam(.bai), reference.fasta(.fai)

    Returns True on success, False on failure (the run as a whole can still
    succeed if at least one sample succeeds).
    """
    key = sample["key"]
    label = sample["label"]
    reads = sample["reads"]

    # Per-sample S3 prefix and local output dir.
    sample_prefix = f"{s3_folder}{key}"               # e.g. "AbCd1234/A_2019_1"
    output_dir = os.path.join(base_output_dir, s3_folder, key)
    print(f"[sample:{key}] reads={reads}")
    print(f"[sample:{key}] output_dir={output_dir}")

    # ---- validate this sample's reads ----
    for fastq_path in reads:
        ok, err = validate_fastq_file(fastq_path)
        if not ok:
            print(f"[sample:{key}] FASTQ check failed: {err}")
            return False
    if not reads:
        print(f"[sample:{key}] no reads")
        return False

    gd_file = os.path.join(output_dir, "output", "output.gd")
    if not os.path.exists(gd_file):
        print(f"[sample:{key}] Running breseq.")
        ok = run_breseq_command(
            os.path.dirname(reads[0]),
            reads,
            output_dir,
            poly,
            referenceFile,
            name=(name or label),
            extra_flags=extra_flags,
            threads=BRESEQ_THREADS_PER_JOB,
        )
        if not ok:
            print(f"[sample:{key}] breseq FAILED")
            return False
    else:
        print(f"[sample:{key}] Existing breseq output found; skipping breseq.")

    if not os.path.exists(output_dir):
        print(f"[sample:{key}] output dir missing after run")
        return False

    # ---- mutation JSON ----
    mutation_file = os.path.join(output_dir, "mutation_predictions.json")
    if not os.path.exists(mutation_file):
        generate_mutation_json(output_dir)

    # ---- coverage + averages ----
    coverage_file = run_samtools_command(output_dir)
    averages_file = os.path.join(output_dir, "averages.csv")
    if coverage_file:
        calculate_coverage_averages(coverage_file, output_dir)

    bam_file = os.path.join(output_dir, "data", "reference.bam")
    bai_file = os.path.join(output_dir, "data", "reference.bam.bai")

    # ---- upload this sample's results to <slug>/<key>/... ----
    print(f"[sample:{key}] Uploading results to s3://{results_bucket}/{sample_prefix}/")
    breseq_html_dir = os.path.join(output_dir, "output")
    if os.path.isdir(breseq_html_dir):
        upload_directory_with_mime(breseq_html_dir, results_bucket, f"{sample_prefix}/output")

    # single downloadable tarball of this sample's output
    tar_path = f"{output_dir.rstrip('/')}.tar.gz"
    os.system(f"tar -czf {tar_path} -C {output_dir.rstrip('/')} .")
    if os.path.exists(tar_path):
        upload_file_to_s3(results_bucket, f"{sample_prefix}/", tar_path)

    for f in (mutation_file, coverage_file, averages_file, bam_file, bai_file):
        if f and os.path.exists(f):
            upload_file_to_s3(results_bucket, f"{sample_prefix}/", f)

    contigsFile = os.path.join(output_dir, "data", "reference.fasta")
    if os.path.exists(contigsFile):
        try:
            subprocess.run(["samtools", "faidx", contigsFile], check=True)
        except Exception as e:
            print(f"[sample:{key}] faidx failed: {e}")
        upload_file_to_s3(results_bucket, f"{sample_prefix}/", contigsFile)
        contigsIndex = contigsFile + ".fai"
        if os.path.exists(contigsIndex):
            upload_file_to_s3(results_bucket, f"{sample_prefix}/", contigsIndex)

    print(f"[sample:{key}] DONE")
    return True


def process_s3_folder(s3_folder, input_bucket, results_bucket, local_base_dir, base_output_dir):
    """
    Process one S3 submission folder end-to-end.

    Each SAMPLE in the submission is run as its OWN breseq job. Inputs are read
    from `input_bucket`; results are written to `results_bucket` under a
    PER-SAMPLE sub-prefix:

        <slug>/<sample_key>/output/index.html
        <slug>/<sample_key>/mutation_predictions.json
        <slug>/<sample_key>/...

    A manifest is written to <slug>/samples.json so the web app's results viewer
    can offer a toggle between samples:

        {"samples": [{"label": "A_2019_1", "key": "A_2019_1"},
                     {"label": "ST_2019_2", "key": "ST_2019_2"}]}
    """
    try:
        print(f"Processing S3 folder: {s3_folder}")

        local_folder = os.path.join(local_base_dir, s3_folder)
        print(f"Local folder path: {local_folder}")

        # Inputs come from the INPUT bucket.
        download_s3_folder(input_bucket, s3_folder, local_folder)

        # Extract form data (no email anymore); now returns per-sample specs.
        referenceFile, poly, read_paths, name, extra_flags, app, sample_specs = \
            extract_form_data(local_folder)

        print(f"Reference: {referenceFile} | poly: {poly} | name: {name!r} | flags: {extra_flags}")
        print(f"Samples: {[s['key'] for s in sample_specs]}")

        if not sample_specs:
            print("[FAILED] No samples / reads found in form-data.txt")
            append_seen_folder(failed_log_file_path, s3_folder)
            return False

        if referenceFile == "None" or not os.path.exists(referenceFile):
            print(f"[FAILED] Reference not available: {referenceFile}")
            append_seen_folder(failed_log_file_path, s3_folder)
            return False

        # ---- run each sample independently ----
        succeeded = []
        for sample in sample_specs:
            ok = _process_one_sample(
                sample, s3_folder, results_bucket, referenceFile, poly,
                name, extra_flags, local_base_dir, base_output_dir,
            )
            if ok:
                succeeded.append(sample)

        if not succeeded:
            print(f"[FAILED] All samples failed for {s3_folder}")
            append_seen_folder(failed_log_file_path, s3_folder)
            return False

        # ---- write the samples manifest the viewer reads ----
        manifest = {
            "slug": s3_folder.rstrip("/"),
            "run_name": name,
            "samples": [{"label": s["label"], "key": s["key"]} for s in succeeded],
        }
        manifest_path = os.path.join(base_output_dir, s3_folder, "samples.json")
        os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
        with open(manifest_path, "w") as mf:
            json.dump(manifest, mf, indent=2)
        # upload to <slug>/samples.json  (inline JSON)
        s3_client.upload_file(
            Filename=manifest_path,
            Bucket=results_bucket,
            Key=f"{s3_folder}samples.json",
            ExtraArgs={"ContentType": "application/json", "ContentDisposition": "inline"},
        )
        print(f"[manifest] Wrote samples.json with {len(succeeded)} sample(s)")

        append_seen_folder(log_file_path, s3_folder)
        print(f"Completed processing for folder: {s3_folder} "
              f"({len(succeeded)}/{len(sample_specs)} samples ok)")
        return True

    except Exception as exc:
        print(f"[FAILED] Unexpected error while processing {s3_folder}: {exc}")
        append_seen_folder(failed_log_file_path, s3_folder)
        print(f"[FAILED] Added to failed folders log: {s3_folder}")
        return False


if __name__ == "__main__":
    input_bucket = INPUT_BUCKET        # midauthorbio-breseq-input  (read inputs)
    results_bucket = RESULTS_BUCKET    # midauthorbio-breseq-results (write results)
    local_base_dir = '/home/ark/MAB/breseq/data'
    folders = list_folders_in_bucket(input_bucket)
    base_output_dir = '/home/ark/MAB/breseq/results'
    app = "breseq"

    print(f"[CONFIG] input bucket  = {input_bucket}")
    print(f"[CONFIG] results bucket = {results_bucket}")

    seen_folders = load_seen_folders(log_file_path)
    failed_folders = load_seen_folders(failed_log_file_path)

    new_folders = []

    for s3_folder in folders:
        if s3_folder in seen_folders:
            continue

        if s3_folder in failed_folders:
            print(f"[SKIP] Previously failed folder: {s3_folder}")
            continue

        # A run is "done" once its per-sample manifest exists in the RESULTS
        # bucket (samples.json is written only after at least one sample
        # completed and uploaded).
        existing_result_key = f"{s3_folder}samples.json"

        if s3_key_exists(results_bucket, existing_result_key):
            print(f"[SKIP] Found existing results for {s3_folder} at s3://{results_bucket}/{existing_result_key}")
            append_seen_folder(log_file_path, s3_folder)
            continue

        new_folders.append(s3_folder)

    # Debugging: Check if folders are retrieved
    print(f"Folders found in input bucket: {folders}")
    if new_folders:
        for folder in new_folders:
            print(f"[NEW FOLDER] {folder}")
    else:
        print("[SCAN] No new folders found")

    if new_folders:
        print(
            f"[PARALLEL] Processing {len(new_folders)} folder(s) "
            f"with up to {MAX_PARALLEL_JOBS} concurrent job(s); "
            f"each breseq job uses {BRESEQ_THREADS_PER_JOB} thread(s)."
        )

        completed = 0
        failed = 0

        with ThreadPoolExecutor(max_workers=MAX_PARALLEL_JOBS) as executor:
            future_to_folder = {
                executor.submit(
                    process_s3_folder,
                    s3_folder,
                    input_bucket,
                    results_bucket,
                    local_base_dir,
                    base_output_dir
                ): s3_folder
                for s3_folder in new_folders
            }

            for future in as_completed(future_to_folder):
                s3_folder = future_to_folder[future]
                try:
                    success = future.result()
                except Exception as exc:
                    print(f"[FAILED] Unhandled worker error for {s3_folder}: {exc}")
                    append_seen_folder(failed_log_file_path, s3_folder)
                    success = False

                if success:
                    completed += 1
                    print(f"[PARALLEL] Completed {s3_folder}")
                else:
                    failed += 1
                    print(f"[PARALLEL] Failed {s3_folder}")

        print(f"[PARALLEL] Completed {completed}; failed {failed}.")

    print("All folders processed.")
