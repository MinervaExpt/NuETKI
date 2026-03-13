import datetime
import subprocess

def extract_table(lines, start_key):
    """
    Extracts table rows following a header like 'MC cut summary:'
    Stops at blank lines or constraint sections.
    """
    in_table = False
    rows = []

    for line in lines:
        if start_key in line:
            in_table = True
            continue

        if in_table:
            if (
                not line.strip()
                or "Signal Constraint" in line
                or "Phase Space Constraint" in line
            ):
                break

            if (
                '|' in line
                and 'Cut Name' not in line
                and '-' not in line
            ):
                rows.append(line)

    return rows


sum_POT = True
sum_mc_table = True
sum_data_table = True

marker = "MacroUtil configuration of this macro runEventLoop"

x = datetime.datetime.now()
outFileName = "CutSummary_" + x.strftime("%b_%d_%Y") + ".txt"

# Read list of cut table file paths ONCE
with open("CutSummary_files.txt") as f:
    file_paths = [line.strip() for line in f if line.strip()]

# -----------------------
# POT accumulators
# -----------------------
n_dataFiles = 0
n_MCFiles = 0
n_truthFiles = 0
data_POT = 0.0
mc_POT = 0.0

# -----------------------
# MC table accumulators
# -----------------------
mc_n_cuts = None
mc_cut_names = []
mc_event_sums = []
mc_eff_weighted = []
mc_purity_weighted = []
mc_rel_eff_weighted = []
mc_rel_all_weighted = []

# -----------------------
# Data table accumulators
# -----------------------
data_n_cuts = None
data_cut_names = []
data_event_sums = []
data_rel_sample_weighted = []

# =======================
# MAIN LOOP OVER FILES
# =======================
for idx, file_path in enumerate(file_paths):
    with open(file_path) as f:
        all_lines = list(f)

    # --- find marker ---
    try:
        start_idx = next(
            i for i, line in enumerate(all_lines) if marker in line
        )
    except StopIteration:
        raise RuntimeError(f"Marker not found in file: {file_path}")

    lines = all_lines[start_idx:]

    # -----------------------
    # POT section
    # -----------------------
    if sum_POT:
        n_dataFiles  += int(lines[1].split()[-1])
        n_MCFiles    += int(lines[2].split()[-1])
        n_truthFiles += int(lines[3].split()[-1])
        data_POT     += float(lines[4].split()[-1])
        mc_POT       += float(lines[5].split()[-1])

    # -----------------------
    # Extract table lines once
    # -----------------------
    mc_table_lines = []
    data_table_lines = []

    if sum_mc_table:
        mc_table_lines = extract_table(lines, "MC cut summary:")

    if sum_data_table:
        data_table_lines = extract_table(lines, "Data cut summary:")

    # -----------------------
    # MC table
    # -----------------------
    if sum_mc_table:
        if mc_n_cuts is None:
            mc_n_cuts = len(mc_table_lines)
            mc_event_sums = [0.0] * mc_n_cuts
            mc_eff_weighted = [0.0] * mc_n_cuts
            mc_purity_weighted = [0.0] * mc_n_cuts
            mc_rel_eff_weighted = [0.0] * mc_n_cuts
            mc_rel_all_weighted = [0.0] * mc_n_cuts

        for i, line in enumerate(mc_table_lines):
            parts = [x.strip() for x in line.split('|') if x.strip()]
            if idx == 0:
                mc_cut_names.append(parts[0])

            events = float(parts[1].replace('e+', 'E+'))
            eff = float(parts[2])
            purity = float(parts[3])
            rel_eff = float(parts[4])
            rel_all = float(parts[5])

            mc_event_sums[i] += events
            mc_eff_weighted[i] += eff * events
            mc_purity_weighted[i] += purity * events
            mc_rel_eff_weighted[i] += rel_eff * events
            mc_rel_all_weighted[i] += rel_all * events

    # -----------------------
    # Data table
    # -----------------------
    if sum_data_table:
        if data_n_cuts is None:
            data_n_cuts = len(data_table_lines)
            data_event_sums = [0.0] * data_n_cuts
            data_rel_sample_weighted = [0.0] * data_n_cuts

        for i, line in enumerate(data_table_lines):
            parts = [x.strip() for x in line.split('|') if x.strip()]
            if idx == 0:
                data_cut_names.append(parts[0])

            events = float(parts[1].replace('e+', 'E+'))
            rel_sample = float(parts[2])

            data_event_sums[i] += events
            data_rel_sample_weighted[i] += rel_sample * events

# =======================
# WRITE OUTPUT
# =======================
with open(outFileName, "w") as f:
    if sum_POT:
        f.write(f"{'Total number of data files is':<35} {n_dataFiles}\n")
        f.write(f"{'Total number of MC files is':<35} {n_MCFiles}\n")
        f.write(f"{'Total number of truth files is':<35} {n_truthFiles}\n")
        f.write(f"{'Total Data POT is':<35} {data_POT}\n")
        f.write(f"{'Total MC POT is':<35} {mc_POT}\n\n")

    if sum_mc_table:
        eff_avg = [mc_eff_weighted[i] / mc_event_sums[i] if mc_event_sums[i] > 0 else 0.0 for i in range(mc_n_cuts)]
        purity_avg = [mc_purity_weighted[i] / mc_event_sums[i] if mc_event_sums[i] > 0 else 0.0 for i in range(mc_n_cuts)]
        rel_eff_avg = [mc_rel_eff_weighted[i] / mc_event_sums[i] if mc_event_sums[i] > 0 else 0.0 for i in range(mc_n_cuts)]
        rel_all_avg = [mc_rel_all_weighted[i] / mc_event_sums[i] if mc_event_sums[i] > 0 else 0.0 for i in range(mc_n_cuts)]

        f.write("MC Cut Summary:\n")
        f.write(f"|{'Cut':>57} | {'Events':>10} | {'% Eff':>7} | {'% Purity':>9} | {'Relative % Eff':>14} | {'Relative % All':>14} |\n")
        f.write("|" + "-" * 127 + "|\n")
        for i in range(mc_n_cuts):
            f.write(
                f"|{mc_cut_names[i]:>57} | {mc_event_sums[i]:10.0f} | "
                f"{eff_avg[i]:7.3f} | {purity_avg[i]:9.3f} | "
                f"{rel_eff_avg[i]:14.3f} | {rel_all_avg[i]:14.3f} |\n"
            )
        f.write("\n")

    if sum_data_table:
        rel_sample_avg = [
            data_rel_sample_weighted[i] / data_event_sums[i]
            if data_event_sums[i] > 0 else 0.0
            for i in range(data_n_cuts)
        ]

        f.write("Data Cut Summary:\n")
        f.write(f"|{'Cut':>57} | {'Events':>10} | {'Relative % Sample Left':>22} |\n")
        f.write("|" + "-" * 96 + "|\n")
        for i in range(data_n_cuts):
            f.write(
                f"|{data_cut_names[i]:>57} | {data_event_sums[i]:10.0f} | "
                f"{rel_sample_avg[i]:22.3f} |\n"
            )

subprocess.run(["cat", outFileName])

