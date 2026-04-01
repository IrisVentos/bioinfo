import sys
from collections import Counter
from spectral_convolution import spectral_convolution
from theoretical_spectrum_peptide import linear_spectrum, cyclic_spectrum

def get_candidate_masses(spectrum, m):
    conv = spectral_convolution(spectrum)
    counts = Counter(x for x in conv if 57 <= x <= 200)
    if not counts:
        return []
    sorted_masses = sorted(counts, key=lambda x: counts[x], reverse=True)
    # take top M with ties  
    cutoff_count = counts[sorted_masses[min(m, len(sorted_masses)) - 1]]
    return [mass for mass in sorted_masses if counts[mass] >= cutoff_count]


def score(peptide, spectrum, cyclic=True):
    theo = Counter(cyclic_spectrum(peptide) if cyclic else linear_spectrum(peptide))
    exp  = Counter(spectrum)
    return sum(min(c, exp.get(m, 0)) for m, c in theo.items())

def trim(leaderboard, spectrum, n):
    if not leaderboard:
        return []
    scored = sorted(leaderboard, key=lambda p: score(p, spectrum, cyclic=False), reverse=True)
    if len(scored) <= n:
        return scored
    cutoff = score(scored[n - 1], spectrum, cyclic=False)
    return [p for p in scored if score(p, spectrum, cyclic=False) >= cutoff]

def convolution_cyclopeptide_sequencing(spectrum, m, n):
    candidate_masses = get_candidate_masses(spectrum, m)
    parent_mass  = max(spectrum)
    leaderboard  = [[]]
    leaders      = []
    leader_score = 0

    while leaderboard:
        leaderboard = [p + [mass] for p in leaderboard for mass in candidate_masses]

        next_board = []
        for pep in leaderboard:
            mass = sum(pep)
            if mass == parent_mass:
                s = score(pep, spectrum, cyclic=True)
                if s > leader_score:
                    leaders, leader_score = [pep], s
                elif s == leader_score:
                    leaders.append(pep)
            elif mass < parent_mass:
                next_board.append(pep)

        leaderboard = trim(next_board, spectrum, n)

    return leaders

if __name__ == "__main__":
    input_file = sys.argv[1] if len(sys.argv) > 1 else "bioinfo_genome_sequencing/datasets/dataset_22.txt"
    lines      = open(input_file).read().splitlines()
    m        = int(lines[0].strip())
    n        = int(lines[1].strip())
    spectrum = list(map(int, lines[2].split()))
    leaders = convolution_cyclopeptide_sequencing(spectrum, m, n)
    print("-".join(map(str, leaders[0])))