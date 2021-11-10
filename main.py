import generateTestCasesAlt
import getMinimizersSyncmers
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick

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

# calculates both the ratio of covered mers compared to an entire sequence, as well as
# the ratio of covered mers per kmer.
# this is done by adding the k-value to the total every time two kmer strings don't intersect.
# if they intersect we add the difference between their indices.
# returns the total divided by the original sequence length as well as
# the total divided by the number of kmers
def ratio_covered_seq(kmers, k_size, seq_length):
    nr_of_covered_chars = k_size

    for i in range(1,len(kmers)):
        diff = kmers[i][1] - kmers[i - 1][1]

        if diff < k_size:
            nr_of_covered_chars += diff
        else:
            nr_of_covered_chars += k_size

    return nr_of_covered_chars / seq_length, nr_of_covered_chars / len(kmers)

def get_list_averages(lst, n):
    res = []
    for i in range(len(lst)):
        sum = 0
        for e in lst[i]:
            sum += e
        res.append(sum / n)

    return res

def ratio_to_percent_lst(lst):
    res_lst = []
    for val in lst:
        res_lst.append(val*100)
    return res_lst


def main():
    # presets for easy modification
    seq_len = 2000
    k = 15
    w = k + 10
    # syncmer algo in getMinimizersSyncmers was modified to take t parameter. If t = -1 it will still
    # use the middle position as previously, otherwise it will use the given positio of t to determine
    # s-mer position
    s = 5
    t = -1

    # From eyeing over Edgar's paper he appeared to get the best results at k = 15 with s = 9 and t = 3.
    # Plugging in these parameters appears to make % of preserved kmers slightly bettter,
    # % of covered sequence a lot better, but number of mers covered per kmer worse
    #s = 9
    #t = 3
    n = 5

    # separate the different parts for ease of reading
    # first we will generate the ratios of conserved kmers

    # results for the ratio will be appended to these lists
    mini_preserved_ratio = []
    sync_preserved_ratio = []

    # one iteration for each error level 1-10%
    for error in range(10):
        # temporary lists to save the results of inner iteration, will be appended to the
        # "main" result list each outer iteration
        mini_preserved_ratio_inner = []
        sync_preserved_ratio_inner = []

        # iterate n times to get n different ratios for the same error level
        for i in range(n):
            seq = generateTestCasesAlt.generate_random_sequence_by_length(seq_len)
            seq_error = generate_error_sequence(seq, round((error+1)/100, 2))

            mini = getMinimizersSyncmers.get_kmer_minimizers(seq,k,w)
            mini_error = getMinimizersSyncmers.get_kmer_minimizers(seq_error,k,w)
            mini_preserved_ratio_inner.append(ratio_preserved_kmers(mini, mini_error))

            sync = getMinimizersSyncmers.get_kmer_syncmers(seq,k,s,t)
            sync_error = getMinimizersSyncmers.get_kmer_syncmers(seq_error, k, s, t)
            sync_preserved_ratio_inner.append(ratio_preserved_kmers(sync, sync_error))

        # append inner list to outer
        mini_preserved_ratio.append(mini_preserved_ratio_inner)
        sync_preserved_ratio.append(sync_preserved_ratio_inner)

    # get the average of the results by adding them together and dividing by n
    mini_preserved_ratio_avg = get_list_averages(mini_preserved_ratio,n)
    sync_preserved_ratio_avg = get_list_averages(sync_preserved_ratio,n)



    # in this section we generate the ratio of covered mers covered by
    # the kmers shared from an original read and an error read
    mini_covered_seq = []
    sync_covered_seq = []

    mini_covered_per_kmer = []
    sync_covered_per_kmer = []

    # one iteration for each error level 1-10%
    for error in range(10):

        # temporary lists to save the results of inner iteration, will be appended to the
        # "main" result list each outer iteration
        mini_covered_seq_inner = []
        sync_covered_seq_inner = []

        mini_covered_per_kmer_inner = []
        sync_covered_per_kmer_inner = []

        # iterate n times to get n different ratios for the same error level
        for i in range(n):
            seq = generateTestCasesAlt.generate_random_sequence_by_length(seq_len)
            seq_error = generate_error_sequence(seq, round((error+1)/100, 2))

            mini = getMinimizersSyncmers.get_kmer_minimizers(seq, k, w)
            mini_error = getMinimizersSyncmers.get_kmer_minimizers(seq_error, k, w)
            mini_shared = get_shared_kmers(mini,mini_error)
            # calculate coverage was modified to return both the ratio of covered sequence
            # as well as ratio of coverage per kmer. Save them temporarily in *_res1, *_res2
            # to append afterwards
            mini_res1, mini_res2 = ratio_covered_seq(mini_shared, k, seq_len)
            mini_covered_seq_inner.append(mini_res1)
            mini_covered_per_kmer_inner.append(mini_res2)

            sync = getMinimizersSyncmers.get_kmer_syncmers(seq, k, s, t)
            sync_error = getMinimizersSyncmers.get_kmer_syncmers(seq_error, k, s, t)
            sync_shared = get_shared_kmers(sync,sync_error)
            sync_res1, sync_res2 = ratio_covered_seq(sync_shared, k, seq_len)
            sync_covered_seq_inner.append(sync_res1)
            sync_covered_per_kmer_inner.append(sync_res2)


        # append inner list to outer
        mini_covered_seq.append(mini_covered_seq_inner)
        sync_covered_seq.append(sync_covered_seq_inner)

        mini_covered_per_kmer.append(mini_covered_per_kmer_inner)
        sync_covered_per_kmer.append(sync_covered_per_kmer_inner)

    # get the average of the results by adding them together and dividing by n
    mini_coverage_avg = get_list_averages(mini_covered_seq,n)
    sync_coverage_avg = get_list_averages(sync_covered_seq, n)

    mini_covered_per_kmer_avg = get_list_averages(mini_covered_per_kmer,n)
    sync_covered_per_kmer_avg = get_list_averages(sync_covered_per_kmer, n)

    x = []
    for i in range(10):
        x.append(str(i+1)+"%")
    xi = list(range(len(x)))

    plt.plot(ratio_to_percent_lst(mini_preserved_ratio_avg), marker='o', color='blue', linestyle='None')
    plt.plot(ratio_to_percent_lst(sync_preserved_ratio_avg), marker='o', color='red', linestyle='None')
    plt.ylim([0,100])
    plt.xticks(xi, x)
    plt.legend(["Minimizers", "Syncmers"])
    plt.title("Comparison of preserved k-mers between minimizer\n and syncmer methods at different error rates")
    plt.xlabel("Error rate")
    plt.ylabel("Percentage of preserved mers")
    plt.savefig("preserved_kmers_mini_vs_sync")
    plt.show()

    plt.plot(ratio_to_percent_lst(mini_coverage_avg), marker='o', color='blue', linestyle='None')
    plt.plot(ratio_to_percent_lst(sync_coverage_avg), marker='o', color='red', linestyle='None')
    plt.ylim([0,100])
    plt.xticks(xi, x)
    plt.legend(["Minimizers", "Syncmers"])
    plt.title(
        "Comparison of percentage of covered sequence between minimizer\n and syncmer methods at different error rates")
    plt.xlabel("Error rate")
    plt.ylabel("Percentage of covered sequence")
    plt.savefig("covered_seq_minis_vs_syncs")
    plt.show()

    plt.plot(mini_covered_per_kmer_avg, marker='o', color='blue', linestyle='None')
    plt.plot(sync_covered_per_kmer_avg, marker='o', color='red', linestyle='None')
    plt.xticks(xi, x)
    plt.legend(["Minimizers", "Syncmers"])
    plt.title(
        "Comparison of percentage of covered sequence PER KMER between minimizer\n and syncmer methods at different error rates")
    plt.xlabel("Error rate")
    plt.ylabel("Number of mers covered per kmer")
    plt.savefig("covered_per_kmer_minis_vs_syncs")
    plt.show()




main()
