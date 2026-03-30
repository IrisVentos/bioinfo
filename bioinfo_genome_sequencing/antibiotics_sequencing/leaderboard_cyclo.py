import sys
from collections import Counter
from peptide_mass_counting import AMINO_ACID_MASS 
from theoretical_spectrum_peptide import linear_spectrum, cyclic_spectrum

MASSES = sorted(set(AMINO_ACID_MASS.values()))

def score(peptide, spectrum, cyclic=True):
    theo = Counter(cyclic_spectrum(peptide) if cyclic else linear_spectrum(peptide))
    exp  = Counter(spectrum)
    return sum(min(c, exp.get(m, 0)) for m, c in theo.items())

def trim(leaderboard, spectrum, n):
    scored = sorted(leaderboard, key=lambda p: score(p, spectrum, cyclic=False), reverse=True)
    if len(scored) <= n:
        return scored
    cutoff = score(scored[n - 1], spectrum, cyclic=False)
    return [p for p in scored if score(p, spectrum, cyclic=False) >= cutoff]

def leaderboard_sequencing(spectrum, n):
    parent_mass = max(spectrum)
    leaderboard = [[]]          # peptides stored as lists of masses
    leader      = []
    leader_score = 0

    while leaderboard:
        # Expand: append every candidate mass to every peptide
        leaderboard = [p + [m] for p in leaderboard for m in MASSES]

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
    input_file = sys.argv[1] if len(sys.argv) > 1 else "bioinfo_genome_sequencing/datasets/dataset_19.txt"
    lines      = open(input_file).read().splitlines()
    n          = int(lines[0].strip())
    spectrum   = list(map(int, lines[1].split()))
    leaders    = leaderboard_sequencing(spectrum, n)
    for pep in leaders:
        print("-".join(map(str, pep)))