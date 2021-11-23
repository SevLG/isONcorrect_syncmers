import sys

import generateTestCasesAlt
import getMinimizersSyncmers
from matplotlib import pyplot as plt
from collections import deque, defaultdict
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


#
#
# pasting functions needed for interval score evaluations, so that I don't need to rename
# a bunch of things
#
#
def get_kmer_minimizers(seq, k_size, w_size):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    window_kmers = deque([seq[i:i + k_size] for i in range(w + 1)])
    curr_min = min(window_kmers)
    minimizers = [(curr_min, list(window_kmers).index(curr_min))]

    for i in range(w + 1, len(seq) - k_size):
        new_kmer = seq[i:i + k_size]
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer:
            curr_min = min(window_kmers)
            minimizers.append((curr_min, list(window_kmers).index(curr_min) + i - w))

        # Previous minimizer still in window, we only need to compare with the recently added kmer
        elif new_kmer < curr_min:
            curr_min = new_kmer
            minimizers.append((curr_min, i))

    return minimizers

def get_kmer_syncmers(seq, k_size, s_size, t = -1):
    w = k_size - s_size

    # get t, the position of s-mer
    # t is chosen to be in the middle of k-mer for chosen syncmer
    if t < 0:
        diff=k_size-s_size
        if diff %2==0:
            t=diff/2
        else:
            t=(diff+1)/2
        t -= 1

    syncmers = []
    # get list of all s-mers in first k-mer
    kmer_smers = deque([seq[i:i + s_size] for i in range(w + 1)])
    for i in range(len(seq) - k_size):
        # add new syncmer to list if its smallest s-mer is at place t
        if list(kmer_smers).index(min(kmer_smers)) == t:
            syncmers.append((seq[i:i+k_size], i))
        # move the window one step to the right by popping the leftmost
        # s-mer and adding one to the right
        kmer_smers.popleft()
        kmer_smers.append(seq[i+k_size-s_size+1:i+k_size+1])

    return syncmers

def get_minimizers_and_positions(reads, w, k, sync = False, s = -1, t = -1):
    # 1. homopolymenr compress read and obtain minimizers
    M = {}
    for r_id in range(len(reads)):
        seq = reads[r_id]
        if not sync:
            minimizers = get_kmer_minimizers(seq, k, w)
        else:
            minimizers = get_kmer_syncmers(seq, k, s, t)

        M[r_id] = minimizers

    return M


from array import array


def get_minimizer_combinations_database(reads, M, k, x_low, x_high):  # generates the Minimizer combinations
    # M2 = defaultdict(lambda: defaultdict(list))
    M2 = defaultdict(lambda: defaultdict(lambda: array("I")))
    tmp_cnt = 0
    forbidden = 'A' * k
    for r_id in M:
        minimizers = M[r_id]
        for (m1, p1), m1_curr_spans in minimizers_comb_iterator(minimizers, k, x_low, x_high):
            for (m2, p2) in m1_curr_spans:
                if m2 == m1 == forbidden:
                    continue

                tmp_cnt += 1
                # t = array('I', [r_id, p1, p2])
                # M2[m1][m2].append( t )
                # M2[m1][m2].append((r_id, p1, p2))


                M2[m1][m2].append(r_id)
                M2[m1][m2].append(p1)
                M2[m1][m2].append(p2)

    #print(tmp_cnt, "MINIMIZER COMBINATIONS GENERATED")
    # import time
    # time.sleep(10)
    # sys.exit()

    avg_bundance = 0
    singleton_minimzer = 0
    cnt = 1
    abundants = []
    for m1 in list(M2.keys()):
        for m2 in list(M2[m1].keys()):
            if len(M2[m1][m2]) > 3:
                avg_bundance += len(M2[m1][m2]) // 3
                cnt += 1
            else:
                del M2[m1][m2]
                singleton_minimzer += 1

            if len(M2[m1][m2]) // 3 > len(reads):
                abundants.append((m1, m2, len(M2[m1][m2]) // 3))
                if m2 == forbidden:  # poly A tail
                    del M2[m1][m2]
    # for m1,m2,ab in sorted(abundants, key=lambda x: x[2], reverse=True):
    # print("Too abundant:", m1, m2, ab, len(reads))

    #print("Average abundance for non-unique minimizer-combs:", avg_bundance / float(cnt))
    #print("Number of singleton minimizer combinations filtered out:", singleton_minimzer)

    return M2


def minimizers_comb_iterator(minimizers, k, x_low, x_high):
    # print("read")
    for i, (m1, p1) in enumerate(minimizers[:-1]):
        m1_curr_spans = []
        for j, (m2, p2) in enumerate(minimizers[i + 1:]):
            if x_low < p2 - p1 and p2 - p1 <= x_high:
                m1_curr_spans.append((m2, p2))
                # yield (m1,p1), (m2, p2)
            elif p2 - p1 > x_high:
                break
        yield (m1, p1), m1_curr_spans[::-1]


def fill_p2(p, all_intervals_sorted_by_finish):
    stop_to_max_j = {stop: j for j, (start, stop, w, _) in enumerate(all_intervals_sorted_by_finish)}
    all_choord_to_max_j = []
    j_max = 0
    for i in range(0, all_intervals_sorted_by_finish[-1][1]):
        if i in stop_to_max_j:
            j_max = stop_to_max_j[i]

        all_choord_to_max_j.append(j_max)

    for j, (start, stop, w, _) in enumerate(all_intervals_sorted_by_finish):
        j_max = all_choord_to_max_j[start-1]
        p.append(j_max)
    return p


def solve_WIS(all_intervals_sorted_by_finish):
    # print("instance size", len(all_intervals_sorted_by_finish))
    # p = [None]
    # fill_p(p, all_intervals_sorted_by_finish)
    p = [None]
    fill_p2(p, all_intervals_sorted_by_finish)
    # if p != p2:
    #     print(p)
    #     print(p2)
    # assert p == p2

    v = [None] + [w * (stop - start) for (start, stop, w, _) in all_intervals_sorted_by_finish]
    OPT = [0]
    for j in range(1, len(all_intervals_sorted_by_finish) + 1):
        OPT.append(max(v[j] + OPT[p[j]], OPT[j - 1]))

    # assert len(p) == len(all_intervals_sorted_by_finish) + 1 == len(v) == len(OPT)

    # Find solution
    opt_indicies = []
    j = len(all_intervals_sorted_by_finish)
    score=0
    while j >= 0:
        if j == 0:
            break
        if v[j] + OPT[p[j]] > OPT[j - 1]:
            opt_indicies.append(
                j - 1)# we have shifted all indices forward by one so we neew to reduce to j -1 because of indexing in python works

            j = p[j]
        else:
            j -= 1
    #OPT should contain the values of the intervals. We take the highest score each read gets
    score = max(OPT)
    #print("OPT",OPT)
    return opt_indicies,score

from itertools import zip_longest

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def add_items(seqs, r_id, p1, p2):
    seqs.append(r_id)
    seqs.append(p1)
    seqs.append(p2)

def find_most_supported_span(r_id, m1, p1, m1_curr_spans, minimizer_combinations_database, reads, all_intervals, k_size,
                             tmp_cnt, read_complexity_cnt, already_computed):
    seq = reads[r_id]
    for (m2, p2) in m1_curr_spans:
        #print(p1,p2)
        #print(m1,m2)
        relevant_reads = minimizer_combinations_database[m1][m2]
        #print(relevant_reads)
        seqs = array("I")  # {} #defaultdict(list)
        added_strings = {}
        locations = {}
        # not_added_strings = set()
        if len(relevant_reads) // 3 >= 3:
            # cnt += 1
            ref_seq = seq[p1: p2 + k_size]
            add_items(seqs, r_id, p1, p2)
            locations[0] = len(seqs) - 3
            added_strings[ref_seq] = 0
            reads_visited = {}
            for relevant_read_id, pos1, pos2 in grouper(relevant_reads, 3):  # relevant_reads:
                if r_id == relevant_read_id:
                    continue

                read_seq = reads[relevant_read_id][1][pos1: pos2 + k_size]
                # read_qual = reads[relevant_read_id][2][pos1: pos2 + k_size]

                if read_seq == ref_seq:
                    # seqs[relevant_read_id] = (pos1, pos2)
                    add_items(seqs, relevant_read_id, pos1, pos2)
                    locations[relevant_read_id] = len(seqs) - 3
                    reads_visited[relevant_read_id] = 0
                    already_computed[relevant_read_id] = (p1, p2, pos1, pos2, 0)
                    continue
                elif relevant_read_id in reads_visited:
                    # print("Prev:", reads_visited[relevant_read_id])
                    # print("Act:", edlib_alignment(ref_seq, read_seq, p_error_sum_thresh*len(ref_seq)) )
                    pass
                # Implement if we see this to recompute all the aligments exact ed here instead!! Thats the only way to guarantee exactly the same
                # or maybe use this traceback to get exact: https://github.com/Martinsos/edlib/pull/132#issuecomment-522258271
                elif read_seq in added_strings:  # == ref_seq:
                    # seqs[relevant_read_id] = (pos1, pos2)
                    add_items(seqs, relevant_read_id, pos1, pos2)
                    locations[relevant_read_id] = len(seqs) - 3
                    reads_visited[relevant_read_id] = added_strings[read_seq]
                    already_computed[relevant_read_id] = (p1, p2, pos1, pos2, added_strings[read_seq])
                    continue

                elif relevant_read_id in already_computed:
                    curr_ref_start, curr_ref_end, curr_read_start, curr_read_end, curr_ed = already_computed[
                        relevant_read_id]

                tmp_cnt += 1

            all_intervals.append((p1 + k_size, p2, len(seqs) // 3, seqs))
    del seqs
    return tmp_cnt, read_complexity_cnt


def get_opt_indicies_and_score(seq, nr_of_err_seq, error, k, w = -1, s = -1, t = -1):
    reads = [seq]
    for i in range(nr_of_err_seq):
        reads.append(generate_error_sequence(seq, error))

    k_size = k
    x_high = 80
    x_low = 14
    if w > 0:
        minimizer_database = get_minimizers_and_positions(reads, w, k_size)
    else:
        minimizer_database = get_minimizers_and_positions(reads, w, k_size, True, s, t)

    minimizer_combinations_database = get_minimizer_combinations_database(reads, minimizer_database, k_size, x_low,
                                                                          x_high)
    tmp_cnt = 0
    opt_indices_total = []
    score_total = []
    all_intervals_for_graph = {}
    for r_id in range(len(reads)):  # , reverse=True):

        read_min_comb = [((m1, p1), m1_curr_spans) for (m1, p1), m1_curr_spans in
                         minimizers_comb_iterator(minimizer_database[r_id], k_size, x_low, x_high)]

        previously_corrected_regions = defaultdict(list)

        seq = reads[r_id]

        if previously_corrected_regions[r_id]:
            read_previously_considered_positions = set(
                [tmp_pos for tmp_p1, tmp_p2, w_tmp, _ in previously_corrected_regions[r_id] for tmp_pos in
                 range(tmp_p1, tmp_p2)])
            group_id = 0
            pos_group = {}
            sorted_corr_pos = sorted(read_previously_considered_positions)
            for p1, p2 in zip(sorted_corr_pos[:-1], sorted_corr_pos[1:]):
                if p2 > p1 + 1:
                    pos_group[p1] = group_id
                    group_id += 1
                    pos_group[p2] = group_id
                else:
                    pos_group[p1] = group_id
            if p2 == p1 + 1:
                pos_group[p2] = group_id
        else:
            read_previously_considered_positions = set()
            pos_group = {}

        already_computed = {}
        read_complexity_cnt = 0

        all_intervals = []
        prev_visited_intervals = []

        for (m1, p1), m1_curr_spans in read_min_comb:
            # If any position is not in range of current corrections: then correct, not just start and stop
            not_prev_corrected_spans = [(m2, p2) for (m2, p2) in m1_curr_spans if not (
                    p1 + k_size in read_previously_considered_positions and p2 - 1 in read_previously_considered_positions)]
            set_not_prev = set(not_prev_corrected_spans)
            not_prev_corrected_spans2 = [(m2, p2) for (m2, p2) in m1_curr_spans if
                                         (m2, p2) not in set_not_prev and (
                                                 p1 + k_size in pos_group and p2 - 1 in pos_group and pos_group[
                                             p1 + k_size] != pos_group[p2 - 1])]
            not_prev_corrected_spans += not_prev_corrected_spans2

            if not_prev_corrected_spans:  # p1 + k_size not in read_previously_considered_positions:
                tmp_cnt, read_complexity_cnt = find_most_supported_span(r_id, m1, p1, not_prev_corrected_spans,
                                                                        minimizer_combinations_database, reads,
                                                                        all_intervals, k_size, tmp_cnt,
                                                                        read_complexity_cnt, already_computed)

        # add prev_visited_intervals to intervals to consider
        all_intervals.extend(prev_visited_intervals)

        if previously_corrected_regions[r_id]:  # add previously corrected regions in to the solver
            all_intervals.extend(previously_corrected_regions[r_id])
            del previously_corrected_regions[r_id]

        if not all_intervals:
            # eprint("Found nothing to correct")
            corrected_seq = seq
        else:
            all_intervals.sort(key=lambda x: x[1])
            # print([www for (_, _,  www, _)  in all_intervals])
            # TODO:these are the infos we want: opt_indices gives the indices of the intervals we chose, score gives the final score for each read
            opt_indicies, score = solve_WIS(
                all_intervals)
            opt_indices_total.append(opt_indicies)
            score_total.append(score)
    return opt_indices_total,score_total
#
#
#
#
#





def main():
    # presets for easy modification
    seq_len = 2000
    k = 15
    w = k + 10
    # syncmer algo in getMinimizersSyncmers was modified to take t parameter. If t = -1 it will still
    # use the middle position as previously, otherwise it will use the given positio of t to determine
    # s-mer position
    s = 9
    t = 3

    # From eyeing over Edgar's paper he appeared to get the best results at k = 15 with s = 9 and t = 3.
    # Plugging in these parameters appears to make % of preserved kmers slightly bettter,
    # % of covered sequence a lot better, but number of mers covered per kmer worse
    #s = 9
    #t = 3
    n = 5

    mini_sync_ratio = []

    # separate the different parts for ease of reading
    # first we will generate the ratios of conserved kmers

    # results for the ratio will be appended to these lists
    mini_preserved_ratio = []
    sync_preserved_ratio = []

    opt_indices_mini = []
    score_mini = []
    opt_indices_sync = []
    score_sync = []



    # one iteration for each error level 1-10%
    for error in range(10):

        error_rate = round((error+1)/100,2)
        # temporary lists to save the results of inner iteration, will be appended to the
        # "main" result list each outer iteration
        mini_preserved_ratio_inner = []
        sync_preserved_ratio_inner = []

        seq = generateTestCasesAlt.generate_random_sequence_by_length(seq_len)

        opt_indices_mini_inner, score_mini_inner = get_opt_indicies_and_score(seq, 20, error_rate, k, w)
        opt_indices_sync_inner, score_sync_inner = get_opt_indicies_and_score(seq, 20, error_rate, k, s=s,t=t)

        opt_indices_mini.append(opt_indices_mini_inner)
        opt_indices_sync.append(opt_indices_sync_inner)
        score_mini.append(score_mini_inner)
        score_sync.append(score_sync_inner)


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


            nr_of_mini = len(mini)
            nr_of_sync = len(sync)
            mini_sync_ratio.append(nr_of_sync / nr_of_mini)

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

    mini_sync_ratio_avg = 0
    for ratio in mini_sync_ratio:
        mini_sync_ratio_avg += ratio

    # the average is taken n times for each error rate, so 10*n times
    mini_sync_ratio_avg = mini_sync_ratio_avg / (10*n)

    print(mini_sync_ratio_avg)

    print(opt_indices_mini[0])
    print(score_mini[0])
    print(opt_indices_sync[0])
    print(score_sync[0])

    x = []
    for i in range(10):
        x.append(str(i+1)+"%")
    xi = list(range(len(x)))

    '''
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

    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(mini_preserved_ratio)
    plt.ylim([0, 1])
    plt.xticks(xi, x)
    plt.title("The spread of preserved minimizers per error rate")
    plt.xlabel("Error rate")
    plt.ylabel("Ratio of preserved kmers")
    plt.savefig("boxplot_mini")
    plt.show()

    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(sync_preserved_ratio)
    plt.ylim([0, 1])
    plt.xticks(xi, x)
    plt.title("The spread of preserved syncmers per error rate")
    plt.xlabel("Error rate")
    plt.ylabel("Ratio of preserved mers")
    plt.savefig("boxplot_sync")
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
        "Comparison of number of mers covered per k-mer between minimizer\n and syncmer methods at different error rates")
    plt.xlabel("Error rate")
    plt.ylabel("Number of mers covered per kmer")
    plt.savefig("covered_per_kmer_minis_vs_syncs")
    plt.show()
    '''





main()
