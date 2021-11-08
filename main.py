import generateTestCasesAlt
import getMinimizersSyncmers
from matplotlib import pyplot as plt

# searches through mer_list and compares a mer with all mers in the list, returns true if found
# otherwise false
def mer_exists(mer, mer_list):
    return (mer in [tup[0] for tup in mer_list])

def get_shared_kmers(s1,s2):
    shared = []
    for i in range(len(s1)):
        if mer_exists(s1[i][0],s2):
            shared.append(s1[i])
    return shared

# generate and return a list of error filled reads with error rates 0.01-0.10 (increments of 0.01)
def generate_error_reads_list(seq):
    isoforms = {}
    acc = "@sim|correct|full"
    isoforms[acc] = seq
    acc = "@sim|correct|error"
    isoforms[acc] = seq
    reads = []
    for error in range(10):
        reads.append(generateTestCasesAlt.simulate_reads(isoforms, round(1-(error+1)/100,2))["@sim|correct|error"][0])
    return reads

# count the number of shared mers between two sets of k-mer tuples (consisting of k-mer sequence and index)
def get_nr_of_shared_mers(s1,s2):
    shared_mers = 0
    for mer in s1:
        if mer_exists(mer[0], s2):
            shared_mers += 1
    return shared_mers

def get_list_of_minimizers(sim_reads, k, w):
    sim_reads_minis = []
    for read in sim_reads:
        sim_reads_minis.append(getMinimizersSyncmers.get_kmer_minimizers(read, k, w))
    return sim_reads_minis

def get_list_of_syncmers(sim_reads, k, s):
    sim_reads_syncs = []
    for read in sim_reads:
        sim_reads_syncs.append(getMinimizersSyncmers.get_kmer_syncmers(read, k, s))
    return sim_reads_syncs

def test_setup(seq_len, k, w, s):
    seq = generateTestCasesAlt.generate_random_sequence_by_length(seq_len)
    mini = getMinimizersSyncmers.get_kmer_minimizers(seq, k, w)
    sync = getMinimizersSyncmers.get_kmer_syncmers(seq, k, s)

    sim_reads = generate_error_reads_list(seq)
    sim_reads_minis = get_list_of_minimizers(sim_reads, k, w)
    sim_reads_syncs = get_list_of_syncmers(sim_reads, k, s)

    return seq, sim_reads, mini, sync, sim_reads_minis, sim_reads_syncs

def ratio_preserved_mers_per_error_rate(kmers, sim_read_kmers):
    nr_shared = []
    for read in sim_read_kmers:
        nr_shared.append(get_nr_of_shared_mers(kmers,read))

    ratios = []
    for i in range(len(nr_shared)):
        ratios.append(nr_shared[i]/((len(kmers)+len(sim_read_kmers[i]))/2))

    return ratios

def lst_sum(lsts):
    res = [0 for i in range(len(lsts[0]))]

    for lst in lsts:
        for i in range(len(lst)):
            res[i] += lst[i]

    return res

def calculate_coverage(k_mers, k_size, seq_length):
    nr_of_covered_chars = k_size
    for i in range(1,len(k_mers)):
        diff = k_mers[i][1] - k_mers[i - 1][1]
        if diff < k_size:
            nr_of_covered_chars += diff
        else:
            nr_of_covered_chars += k_size

    return nr_of_covered_chars / seq_length

def main():
    seq_len = 2000
    k = 15
    w = k + 10
    s = 5
    seq, sim_reads, mini, sync, sim_reads_minis, sim_reads_syncs = test_setup(seq_len, k, w, s)
    for read in sim_reads_syncs:
        print(read)
    '''
    seq_len = 2000
    k = 15
    w = k+10
    s = 5
    n = 5
    ratios_mini = []
    ratios_sync = []
    avg_cover_mini = []
    avg_cover_sync = []
    for i in range(n):
        seq, sim_reads, mini, sync, sim_reads_minis, sim_reads_syncs = test_setup(seq_len, k, w, s)
        ratios_mini.append(ratio_preserved_mers_per_error_rate(mini, sim_reads_minis))
        ratios_sync.append(ratio_preserved_mers_per_error_rate(sync, sim_reads_syncs))

        shared_minis = []
        shared_syncs = []
        for j in range(len(sim_reads_minis)):
            shared_minis.append(get_shared_kmers(sim_reads_minis[j], mini))
            shared_syncs.append(get_shared_kmers(sim_reads_syncs[j], sync))
        cover_mini = []
        cover_sync = []
        for l in range(len(shared_minis)):
            cover_mini.append(calculate_coverage(shared_minis[l], k, seq_len))
            cover_sync.append(calculate_coverage(shared_syncs[l], k, seq_len))
        avg_cover_mini.append(cover_mini)
        avg_cover_sync.append(cover_sync)

    mini_total = lst_sum(ratios_mini)
    sync_total = lst_sum(ratios_sync)
    mini_total_coverage = lst_sum(avg_cover_mini)
    sync_total_coverage = lst_sum(avg_cover_sync)

    for i in range(len(mini_total)):
        mini_total[i] = mini_total[i] / n
        sync_total[i] = sync_total[i] / n
        mini_total_coverage[i] = mini_total_coverage[i] / n
        sync_total_coverage[i] = sync_total_coverage[i] / n

    x = []
    for i in range(10):
        x.append(round(1 - (i + 1) / 100, 2))
    xi = list(range(len(x)))
    plt.plot(mini_total, marker='o', color='blue', linestyle='None')
    plt.plot(sync_total, marker='o', color='red', linestyle='None')
    # plt.hlines(y=0, xmin=-0.1, xmax=1, color=colors[4], linestyles=(0, (1, 5)))
    # plt.xlim([-0.1, 1])
    # plt.ylim([min(errors) - 0.1, max(errors) + 0.1])
    plt.xticks(xi, x)
    plt.legend(["Minimizers", "Syncmers"])
    plt.title("Comparison of preserved k-mers between minimizer\n and syncmer methods at different error rates")
    plt.savefig("preserved_kmers_mini_vs_sync")
    plt.show()

    plt.plot(mini_total_coverage, marker='o', color='blue', linestyle='None')
    plt.plot(sync_total_coverage, marker='o', color='red', linestyle='None')
    # plt.hlines(y=0, xmin=-0.1, xmax=1, color=colors[4], linestyles=(0, (1, 5)))
    # plt.xlim([-0.1, 1])
    # plt.ylim([min(errors) - 0.1, max(errors) + 0.1])
    plt.xticks(xi, x)
    plt.legend(["Minimizers", "Syncmers"])
    plt.title(
        "Comparison of percentage of covered sequence between minimizer\n and syncmer methods at different error rates")
    plt.savefig("covered_seq_minis_vs_syncs")
    plt.show()
    '''

main()

