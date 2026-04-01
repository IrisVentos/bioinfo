import sys
from collections import Counter
from theoretical_spectrum_peptide import linear_spectrum, cyclic_spectrum

def extended_mass_table():
    """Returns a dict mapping char -> int mass for ASCII 57..200."""
    return {chr(i): i for i in range(57, 201)}

# Pre-compute the mass list once (integers)
EXTENDED_MASSES = list(extended_mass_table().values())  # [57, 58, ..., 200]

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
    parent_mass  = max(spectrum)
    leaderboard  = [[]]          # each peptide is a list of int masses
    leaders      = []
    leader_score = 0

    while leaderboard:
        # Expand: add each candidate mass to every current peptide
        leaderboard = [p + [m] for p in leaderboard for m in EXTENDED_MASSES]

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
    input_file = sys.argv[1] if len(sys.argv) > 1 else "bioinfo_genome_sequencing/datasets/dataset_20.txt"
    lines      = open(input_file).read().splitlines()
    n          = int(lines[0].strip())
    spectrum   = list(map(int, lines[1].split()))
    leaders    = leaderboard_sequencing(spectrum, n)
    for pep in leaders:
        print("-".join(map(str, pep)))