from collections import deque

def get_kmer_minimizers(seq, k_size, w_size):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    window_kmers = deque([seq[i:i+k_size] for i in range(w +1)])
    curr_min = min(window_kmers)
    minimizers = [ (curr_min, list(window_kmers).index(curr_min)) ]

    for i in range(w+1,len(seq) - k_size):
        new_kmer = seq[i:i+k_size]
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer:
            curr_min = min(window_kmers)
            minimizers.append( (curr_min, list(window_kmers).index(curr_min) + i - w ) )

        # Previous minimizer still in window, we only need to compare with the recently added kmer
        elif new_kmer < curr_min:
            curr_min = new_kmer
            minimizers.append( (curr_min, i) )

    return minimizers

def get_kmer_syncmers_window_exclusive(seq, k_size, s_size):
    w = k_size - s_size

    # get t, the position of s-mer
    # t is chosen to be in the middle of k-mer for chosen syncmer
    diff=k_size-s_size
    if diff %2==0:
        t=diff/2
    else:
        t=(diff+1)/2
    t -= 1
    syncmers = []
    # get list of all s-mers in first k-mer
    kmer_smers = deque([seq[i:i + s_size] for i in range(w + 1)])
    # we only want the first found syncmer, so we use syncmer_found to see if
    # one has been found in a given window
    syncmer_found = False
    for i in range(len(seq) - k_size):
        # keeping track of window index, 0 == new window
        window_index = i % w
        if window_index == 0:
            syncmer_found = False
        # add new syncmer to list if its smallest s-mer is at place t and
        # another syncmer has not been added in this window
        if list(kmer_smers).index(min(kmer_smers)) == t and not syncmer_found:
            syncmers.append((seq[i:i+k_size], i))
            syncmer_found = True
        # move the window one step to the right by popping the leftmost
        # s-mer and adding one to the right
        kmer_smers.popleft()
        kmer_smers.append(seq[i+k_size-s_size+1:i+k_size+1])

    return syncmers

def get_kmer_syncmers(seq, k_size, s_size):
    w = k_size - s_size

    # get t, the position of s-mer
    # t is chosen to be in the middle of k-mer for chosen syncmer
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




def minimizer_syncmer_comparison():
    k_size = 9
    w_size = k_size + 2
    # seq = "CTGACCGTAC"
    # seq = "GTATCGGCATTACTGACACGAATCGTCAG"
    seq = "CTGGAATCTTAACTGACTGAATCGCTTCGATATATGACTAGCATTCATGC"
    print("Sequence: " + seq)
    minimizers = get_kmer_minimizers(seq,k_size,w_size)
    print("Minimizers:")
    print(minimizers)
    s_size = 4
    syncmers = get_kmer_syncmers(seq,k_size,s_size)
    print("Syncmers:")
    print(syncmers)

