import datetime
import subprocess

sum_mc_table = True
sum_data_table = True
sum_POT = True

#x = datetime.datetime.now()
x = datetime.datetime(2025, 6, 5)
outFileName = "CutSummary_" + x.strftime("%b_%d_%Y") + ".txt"

if sum_POT:
    with open("cutTables.txt") as f:
        file_paths = [line.strip() for line in f if line.strip()]

    n_dataFiles = 0
    n_MCFiles = 0
    n_truthFiles = 0
    data_POT = 0
    mc_POT = 0

    for idx, file_path in enumerate(file_paths):
        with open(file_path) as f:
            lines = [line for line in f]

            n_dataFiles += int(lines[1].strip().split()[-1])
            n_MCFiles += int(lines[2].strip().split()[-1])
            n_truthFiles += int(lines[3].strip().split()[-1])
            data_POT += float(lines[4].strip().split()[-1])
            mc_POT += float(lines[5].strip().split()[-1])

    with open(outFileName, "w") as f:
        # write totals
        f.write(f"{'Total number of data files is':<35} {n_dataFiles}\n")
        f.write(f"{'Total number of MC files is':<35} {n_MCFiles}\n")
        f.write(f"{'Total number of truth files is':<35} {n_truthFiles}\n")
        f.write(f"{'Total Data POT is':<35} {data_POT}\n")
        f.write(f"{'Total MC POT is':<35} {mc_POT}\n\n")
 

if sum_mc_table:
    with open("cutTables.txt") as f:
        file_paths = [line.strip() for line in f if line.strip()]

    n_cuts = None
    cut_names = []
    event_sums = []
    eff_weighted_sums = []
    purity_weighted_sums = []
    relative_eff_weighted_sums = []
    relative_all_weighted_sums = []
    for idx, file_path in enumerate(file_paths):
        with open(file_path) as f:
            # Only keep table lines (contain '|' but not headers or dividers)
            lines = [line for line in f if '|' in line and 'Cut Name' not in line and '-' not in line]

            if n_cuts is None:
                n_cuts = len(lines)
                event_sums = [0.0] * n_cuts
                eff_weighted_sums = [0.0] * n_cuts
                purity_weighted_sums = [0.0] * n_cuts
                relative_eff_weighted_sums = [0.0] * n_cuts
                relative_all_weighted_sums = [0.0] * n_cuts

            for i, line in enumerate(lines):
                parts = [x.strip() for x in line.strip().split('|') if x.strip()]
                if idx == 0:
                    cut_names.append(parts[0])
                events = float(parts[1].replace('e+', 'E+'))
                eff = float(parts[2])
                purity = float(parts[3])
                relative_eff = float(parts[4])
                relative_all = float(parts[5])

                event_sums[i] += events
                eff_weighted_sums[i] += eff * events
                purity_weighted_sums[i] += purity * events
                relative_eff_weighted_sums[i] += relative_eff * events
                relative_all_weighted_sums[i] += relative_all * events

    # Calculate weighted averages
    eff_avg = [eff_weighted_sums[i] / event_sums[i] if event_sums[i] > 0 else 0.0 for i in range(n_cuts)]
    purity_avg = [purity_weighted_sums[i] / event_sums[i] if event_sums[i] > 0 else 0.0 for i in range(n_cuts)]
    relative_eff_avg = [relative_eff_weighted_sums[i] / event_sums[i] if event_sums[i] > 0 else 0.0 for i in range(n_cuts)]
    relative_all_avg = [relative_all_weighted_sums[i] / event_sums[i] if event_sums[i] > 0 else 0.0 for i in range(n_cuts)]

    with open(outFileName, "a") as f:
        f.write("MC Cut Summary: \n")
        # Write combined MC cut table
        f.write(f"|{'Cut':>57} | {'Events':>10} | {'% Eff':>7} | {'% Purity':>9} | {'Relative % Eff':>14} | {'Relative % All':>14} |\n")
        f.write("|" + '-' * 127 + "|\n")
        for i in range(n_cuts):
            f.write(f"|{cut_names[i]:>57} | {event_sums[i]:10.0f} | {eff_avg[i]:7.3f} | {purity_avg[i]:9.3f} | {relative_eff_avg[i]:14.3f} | {relative_all_avg[i]:14.3f} |\n")

        f.write("\n")
        
if sum_data_table:
    # Read list of cut table file paths
    with open("cutTables.txt") as f:
        file_paths = [line.strip() for line in f if line.strip()]

    n_cuts = None
    cut_names = []
    event_sums = []
    relativeSample_weighted_sums = []

    for idx, file_path in enumerate(file_paths):
        with open(file_path) as f:
            # Only keep table lines (contain '|' but not headers or dividers)
            lines = [line for line in f if '|' in line and 'Cut Name' not in line and '-' not in line]
            
            if n_cuts is None:
                n_cuts = len(lines)
                event_sums = [0.0] * n_cuts
                relativeSample_weighted_sums = [0.0] * n_cuts
                
            for i, line in enumerate(lines):
                parts = [x.strip() for x in line.strip().split('|') if x.strip()]
                if idx == 0:
                    cut_names.append(parts[0])
                events = float(parts[1].replace('e+', 'E+'))
                relativeSample = float(parts[2])
                event_sums[i] += events
                relativeSample_weighted_sums[i] += relativeSample * events

    # Calculate weighted averages
    relativeSample_avg = [relativeSample_weighted_sums[i] / event_sums[i] if event_sums[i] > 0 else 0.0 for i in range(n_cuts)]

    with open(outFileName, "a") as f:
        f.write("Data Cut Summary: \n")
        # write combined cut table
        f.write(f"|{'Cut':>57} | {'Events':>10} | {'Relative % Sample Left':>22} |\n")
        f.write('|' + '-' * 96 + "|\n")
        for i in range(n_cuts):
            f.write(f"|{cut_names[i]:>57} | {event_sums[i]:10.0f} | {relativeSample_avg[i]:22.3f} |\n")


subprocess.run(["cat", outFileName])
