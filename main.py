import sys

import generateTestCasesAlt
import getMinimizersSyncmers
from matplotlib import pyplot as plt
import numpy as np; np.random.seed(1)
from intervals_test import get_opt_indicies_and_score

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

# code taken from https://matplotlib.org/2.0.2/examples/pylab_examples/barchart_demo.html
# used to generate barplots side by side from two different datasets
def get_barplot_two_datasets(data1, data2, n_groups):
    index = np.arange(n_groups)
    fig, ax = plt.subplots()
    bar_width = 0.35
    opacity = 0.4
    error_config = {'ecolor': '0.3'}
    rects1 = plt.bar(index, data1, bar_width,
                     alpha=opacity,
                     color='b',
                     error_kw=error_config,
                     label='Minimizers')

    rects2 = plt.bar(index + bar_width, data2, bar_width,
                     alpha=opacity,
                     color='r',
                     error_kw=error_config,
                     label='Syncmers')
    plt.legend()

# used to generate boxplots side by side from two different datasets
def get_boxplot_two_datasets(data1,data2):
    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    pos = [i - 0.2 for i in range(10)]
    bp = ax.boxplot(data1, positions=pos, widths=0.3, patch_artist=True)
    for element in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp[element], color="blue")
    for patch in bp['boxes']:
        patch.set(facecolor="white")
    pos = [i + 0.2 for i in range(10)]
    bp = ax.boxplot(data2, positions=pos, widths=0.3, patch_artist=True)
    for element in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp[element], color="red")
    for patch in bp['boxes']:
        patch.set(facecolor="white")
    leg = plt.legend(["Minimizers", "Syncmers"])
    leg.legendHandles[0].set_color('blue')
    leg.legendHandles[1].set_color('red')


def main():
    # presets for easy modification
    seq_len = 10000
    k = 15
    w = k + 10
    # syncmer algo in getMinimizersSyncmers was modified to take t parameter. If t = -1 it will still
    # use the middle position as previously, otherwise it will use the given positio of t to determine
    # s-mer position
    s = 11
    t = 1
    n = 10

    # separate test for ratio of number of syncmers to number of minimizers
    # should be function
    '''
    nr_of_mini = 0
    nr_of_sync = 0
    for i in range(n):
        seq = generateTestCasesAlt.generate_random_sequence_by_length(seq_len)
        mini = getMinimizersSyncmers.get_kmer_minimizers(seq, k, w)
        sync = getMinimizersSyncmers.get_kmer_syncmers(seq, k, s, t)
        nr_of_mini += len(mini)
        nr_of_sync += len(sync)
    nr_of_mini = nr_of_mini / n
    nr_of_sync = nr_of_sync / n

    print(nr_of_mini)
    print(nr_of_sync)
    print(1-nr_of_sync / nr_of_mini)
    exit()
    '''

    # keep track of number of minimizers and syncmers
    nr_of_mini = 0
    nr_of_sync = 0

    # separate the different parts for ease of reading
    # first we will generate the ratios of conserved kmers

    # results for the ratio will be appended to these lists
    mini_preserved_ratio = []
    sync_preserved_ratio = []

    # not actually used in this analysis
    opt_indices_mini = []
    opt_indices_sync = []

    # results of the isONcorrect interval score, which is used to statistically error correct sequences
    score_mini_avg = []
    score_sync_avg = []

    score_mini_lst = []
    score_sync_lst = []


    # one iteration for each error level 1-10%
    for error in range(10):

        error_rate = round((error+1)/100,2)
        # temporary lists to save the results of inner iteration, will be appended to the
        # "main" result list each outer iteration
        mini_preserved_ratio_inner = []
        sync_preserved_ratio_inner = []

        score_mini = []
        score_sync = []

        nr_of_mini_inner = 0
        nr_of_sync_inner = 0


        # iterate n times to get n different ratios for the same error level
        # these will be averaged to get a somewhat "fair" result statistically
        for i in range(n):
            seq = generateTestCasesAlt.generate_random_sequence_by_length(seq_len)
            seq_error = generate_error_sequence(seq, round((error+1)/100, 2))


            opt_indices_mini_inner, score_mini_inner = get_opt_indicies_and_score(seq, 50, error_rate, k, w)
            opt_indices_sync_inner, score_sync_inner = get_opt_indicies_and_score(seq, 50, error_rate, k, s=s, t=t)
            score_mini.append(sum(score_mini_inner[1:])/50)
            score_sync.append(sum(score_sync_inner[1:])/50)


            mini = getMinimizersSyncmers.get_kmer_minimizers(seq,k,w)
            mini_error = getMinimizersSyncmers.get_kmer_minimizers(seq_error,k,w)
            mini_preserved_ratio_inner.append(ratio_preserved_kmers(mini, mini_error))

            sync = getMinimizersSyncmers.get_kmer_syncmers(seq,k,s,t)
            sync_error = getMinimizersSyncmers.get_kmer_syncmers(seq_error, k, s, t)
            sync_preserved_ratio_inner.append(ratio_preserved_kmers(sync, sync_error))

            nr_of_mini_inner += len(mini)
            nr_of_sync_inner += len(sync)

        # append inner list to outer
        mini_preserved_ratio.append(mini_preserved_ratio_inner)
        sync_preserved_ratio.append(sync_preserved_ratio_inner)

        score_mini_lst.append(score_mini)
        score_sync_lst.append(score_sync)

        score_mini_avg.append(sum(score_mini)/n)
        score_sync_avg.append(sum(score_sync)/n)

        nr_of_mini += nr_of_mini_inner / n
        nr_of_sync += nr_of_sync_inner / n


    # get the average of the results by adding them together and dividing by n
    mini_preserved_ratio_avg = get_list_averages(mini_preserved_ratio,n)
    sync_preserved_ratio_avg = get_list_averages(sync_preserved_ratio,n)

    nr_of_mini = nr_of_mini / 10
    nr_of_sync = nr_of_sync / 10

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

    # print out parameters and avg k-mers which are included in the result
    print("#BP: " + str(seq_len))
    print("k: " + str(k))
    print("w: " + str(w))
    print("s: " + str(s))
    print("t: " + str(t))
    print("n: " + str(n))
    print("avg #minimizers: " + str(round(nr_of_mini,2)))
    print("avg #syncmers: " + str(round(nr_of_sync,2)))
    print("ratio #syncmers/#minimizers: " + str(round(nr_of_sync/nr_of_mini,2)))


    # doing some magic to make the plots work
    x = []
    for i in range(10):
        x.append(str(i+1)+"%")
    xi = list(range(len(x)))

    n_groups = len(mini_preserved_ratio_avg)

    get_barplot_two_datasets(ratio_to_percent_lst(score_mini_avg),ratio_to_percent_lst(score_sync_avg),n_groups)
    plt.xticks(xi, x)
    frame1 = plt.gca()
    frame1.axes.get_yaxis().set_ticks([])
    plt.xlabel("Error rate")
    plt.ylabel("Average score")
    plt.tight_layout()
    plt.savefig("interval_score_mini_vs_sync_bar")
    plt.show()
    
    get_boxplot_two_datasets(score_mini_lst,score_sync_lst)
    plt.xticks(xi, x)
    frame1 = plt.gca()
    frame1.axes.get_yaxis().set_ticks([])
    plt.xlabel("Error rate")
    plt.ylabel("Average score")
    plt.tight_layout()
    plt.savefig("interval_score_mini_vs_sync_box")
    plt.show()

    get_barplot_two_datasets(ratio_to_percent_lst(mini_preserved_ratio_avg),ratio_to_percent_lst(sync_preserved_ratio_avg),n_groups)
    plt.ylim([0,100])
    plt.xticks(xi, x)
    plt.xlabel("Error rate")
    plt.ylabel("Percentage of preserved mers")
    plt.tight_layout()
    plt.savefig("preserved_kmers_mini_vs_sync_bar")
    plt.show()


    get_boxplot_two_datasets(mini_preserved_ratio,sync_preserved_ratio)
    plt.xticks(xi, x)
    plt.ylim([0, 1])
    plt.xlabel("Error rate")
    plt.ylabel("Ratio of preserved k-mers")
    plt.tight_layout()
    plt.savefig("preserved_kmers_mini_vs_sync_box")
    plt.show()

    get_barplot_two_datasets(ratio_to_percent_lst(mini_coverage_avg),ratio_to_percent_lst(sync_coverage_avg),n_groups)
    plt.ylim([0, 100])
    plt.xticks(xi, x)
    plt.xlabel("Error rate")
    plt.ylabel("Covered sequence")
    plt.tight_layout()
    plt.savefig("covered_seq_minis_vs_syncs_bar")
    plt.show()

    get_boxplot_two_datasets(mini_covered_seq,sync_covered_seq)
    plt.ylim([0, 1])
    plt.xticks(xi, x)
    plt.xlabel("Error rate")
    plt.ylabel("Covered sequence ratio")
    plt.tight_layout()
    plt.savefig("covered_seq_minis_vs_syncs_box")
    plt.show()


    get_barplot_two_datasets(mini_covered_per_kmer_avg,sync_covered_per_kmer_avg,n_groups)
    plt.xticks(xi, x)
    frame1 = plt.gca()
    frame1.axes.get_yaxis().set_ticks([])
    plt.xlabel("Error rate")
    plt.ylabel("Number of mers covered per kmer")
    plt.tight_layout()
    plt.savefig("covered_per_kmer_minis_vs_syncs_bar")
    plt.show()

    get_boxplot_two_datasets(mini_covered_per_kmer,sync_covered_per_kmer)
    plt.xticks(xi, x)
    frame1 = plt.gca()
    frame1.axes.get_yaxis().set_ticks([])
    plt.xlabel("Error rate")
    plt.ylabel("Number of basepairs covered per k-mer")
    plt.tight_layout()
    plt.savefig("covered_per_kmer_minis_vs_syncs_box")
    plt.show()






main()
