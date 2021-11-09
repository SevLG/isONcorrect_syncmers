import generateTestCasesAlt
import getMinimizersSyncmers
from matplotlib import pyplot as plt

# searches through kmer_list and compares a k-mer with all k-mers in the list, returns true if found
# otherwise false.
# Used instead of sets since we need to keep all information in the tuples
def kmer_exists(kmer, kmer_list):
    return (kmer in [tup[0] for tup in kmer_list])

# uses kmer_exists to find all kmers in list s1 that is also in list s2.
# only the sequences are compared, index is disregarded.
# only the indices of the first kmer list are saved.
def get_shared_kmers(s1,s2):
    shared = []
    for i in range(len(s1)):
        # kmer_exists takes a string and a list of tuples, so we call s1[i][0] to get only the kmer string
        if kmer_exists(s1[i][0],s2):
            shared.append(s1[i])
    return shared

# count the number of shared kmers using get_shared_kmers function
def get_nr_of_shared_kmers(s1,s2):
    shared_kmers = get_shared_kmers(s1,s2)
    return len(shared_kmers)

# takes a sequence (string) and an error rate between 0.0 and 1.0.
# function creates a dictionary and fills two entries with the sequence.
# it then mutates the sequence located at "@sim|correct|error" according to an error rate
def generate_error_sequence(seq, error):
    # copied from "main" function in generateTestCases, redundant to make a dict in this
    # case but was easiest way to make it work quickly
    isoforms = {}
    acc = "@sim|correct|full"
    isoforms[acc] = seq
    acc = "@sim|correct|error"
    isoforms[acc] = seq
    # 1 - error because simulate_reads interprets 0.99 to mean 1% error
    error_read = generateTestCasesAlt.simulate_reads(isoforms,1-error)["@sim|correct|error"][0]
    return error_read

# calculates the ratio of preserved kmers between a list of kmer tuples and another list of kmer tuples.
# ratio is calculated by dividing number of shared kmers by the average amount of kmers between the
# two kmer lists
def ratio_preserved_kmers(kmers, error_read_kmers):
    nr_shared = get_nr_of_shared_kmers(kmers, error_read_kmers)

    ratio = nr_shared/((len(kmers)+len(error_read_kmers))/2)

    return ratio

# calculates the amount of characters covered by a set of kmers.
# this is done by adding the k-value to the total every time two kmer strings don't intersect.
# if they intersect we add the difference between their indices.
# the total is then divided by the original sequence length
def calculate_coverage(kmers, k_size, seq_length):
    nr_of_covered_chars = k_size

    for i in range(1,len(kmers)):
        diff = kmers[i][1] - kmers[i - 1][1]

        if diff < k_size:
            nr_of_covered_chars += diff
        else:
            nr_of_covered_chars += k_size

    return nr_of_covered_chars / seq_length


def main():
    # presets for easy modification
    seq_len = 5000
    k = 15
    w = k + 10
    s = 5
    # t = 3
    n = 5

    # separate the different parts for ease of reading
    # first we will generate the ratios of conserved kmers

    # results for the ratio will be appended to these lists
    ratio_preserved_minis = []
    ratio_preserved_syncs = []

    # one iteration for each error level 1-10%
    for error in range(10):
        # temporary lists to save the results of inner iteration, will be appended to the
        # "main" result list each outer iteration
        ratio_preserved_minis_inner_list = []
        ratio_preserved_syncs_inner_list = []

        # iterate n times to get n different ratios for the same error level
        for i in range(n):
            seq = generateTestCasesAlt.generate_random_sequence_by_length(seq_len)
            seq_error = generate_error_sequence(seq, round((error+1)/100, 2))

            mini = getMinimizersSyncmers.get_kmer_minimizers(seq,k,w)
            mini_error = getMinimizersSyncmers.get_kmer_minimizers(seq_error,k,w)
            ratio_preserved_minis_inner_list.append(ratio_preserved_kmers(mini, mini_error))

            sync = getMinimizersSyncmers.get_kmer_syncmers(seq,k,s,t)
            sync_error = getMinimizersSyncmers.get_kmer_syncmers(seq_error, k, s, t)
            ratio_preserved_syncs_inner_list.append(ratio_preserved_kmers(sync, sync_error))

        # append inner list to outer
        ratio_preserved_minis.append(ratio_preserved_minis_inner_list)
        ratio_preserved_syncs.append(ratio_preserved_syncs_inner_list)

    # get the average of the results by adding them together and dividing by n
    ratio_preserved_minis_avg = []
    ratio_preserved_syncs_avg = []
    for i in range(10):
        sum = 0
        for e in ratio_preserved_minis[i]:
            sum += e
        ratio_preserved_minis_avg.append(sum/n)

        sum = 0
        for e in ratio_preserved_syncs[i]:
            sum += e
        ratio_preserved_syncs_avg.append(sum / n)



    # in this section we generate the ratio of covered mers covered by
    # the kmers shared from an original read and an error read
    mini_coverage = []
    sync_coverage = []

    # one iteration for each error level 1-10%
    for error in range(10):

        # temporary lists to save the results of inner iteration, will be appended to the
        # "main" result list each outer iteration
        mini_covered_inner_list = []
        sync_covered_inner_list = []
        
        # iterate n times to get n different ratios for the same error level
        for i in range(n):
            seq = generateTestCasesAlt.generate_random_sequence_by_length(seq_len)
            seq_error = generate_error_sequence(seq, round((error+1)/100, 2))

            mini = getMinimizersSyncmers.get_kmer_minimizers(seq, k, w)
            mini_error = getMinimizersSyncmers.get_kmer_minimizers(seq_error, k, w)
            mini_shared = get_shared_kmers(mini,mini_error)
            mini_covered_inner_list.append(calculate_coverage(mini_shared,k,seq_len))

            sync = getMinimizersSyncmers.get_kmer_syncmers(seq, k, s, t)
            sync_error = getMinimizersSyncmers.get_kmer_syncmers(seq_error, k, s, t)
            sync_shared = get_shared_kmers(sync,sync_error)
            sync_covered_inner_list.append(calculate_coverage(sync_shared,k,seq_len))

        # append inner list to outer
        mini_coverage.append(mini_covered_inner_list)
        sync_coverage.append(sync_covered_inner_list)

    # get the average of the results by adding them together and dividing by n
    mini_coverage_avg = []
    sync_coverage_avg = []
    for i in range(10):
        sum = 0
        for e in mini_coverage[i]:
            sum += e
        mini_coverage_avg.append(sum/n)

        sum = 0
        for e in sync_coverage[i]:
            sum += e
        sync_coverage_avg.append(sum / n)


    x = []
    for i in range(10):
        x.append(round(1 - (i + 1) / 100, 2))
    xi = list(range(len(x)))
    plt.plot(ratio_preserved_minis_avg, marker='o', color='blue', linestyle='None')
    plt.plot(ratio_preserved_syncs_avg, marker='o', color='red', linestyle='None')
    # plt.hlines(y=0, xmin=-0.1, xmax=1, color=colors[4], linestyles=(0, (1, 5)))
    # plt.xlim([-0.1, 1])
    plt.ylim([0,1])
    plt.xticks(xi, x)
    plt.legend(["Minimizers", "Syncmers"])
    plt.title("Comparison of preserved k-mers between minimizer\n and syncmer methods at different error rates")
    plt.savefig("preserved_kmers_mini_vs_sync")
    plt.show()

    plt.plot(mini_coverage_avg, marker='o', color='blue', linestyle='None')
    plt.plot(sync_coverage_avg, marker='o', color='red', linestyle='None')
    # plt.hlines(y=0, xmin=-0.1, xmax=1, color=colors[4], linestyles=(0, (1, 5)))
    # plt.xlim([-0.1, 1])
    plt.ylim([0,1])
    plt.xticks(xi, x)
    plt.legend(["Minimizers", "Syncmers"])
    plt.title(
        "Comparison of percentage of covered sequence between minimizer\n and syncmer methods at different error rates")
    plt.savefig("covered_seq_minis_vs_syncs")
    plt.show()




main()
