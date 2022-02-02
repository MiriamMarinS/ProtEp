def occurrences(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count


def occurrences_mismatch(string, sub, mismatch_allowed, epitope_sequences):
    epitopes = []
    for k, v in epitope_sequences.items():
        if v[0] not in epitopes:
            epitopes.append(v[0])

    count = 0
    for i in range(0, len(string) - len(sub) + 1):
        mismatch = 0
        amplicon = string[i:i+len(sub)]

        #if any(epitope == amplicon for epitope in epitopes):
        #    continue
        for ii in range(0, len(sub)):
            if amplicon[ii] != sub[ii]:
                mismatch += 1
        if mismatch == mismatch_allowed:
            count += 1
    return(count)