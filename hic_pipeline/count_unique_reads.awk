# From https://github.com/aidenlab/juicer/blob/e142fbae5a21b152851f65bab11ed051500e3e4b/CPU/juicer.sh#L562
NR == FNR {
    if ($1 ~ /Alignable/) {
        split($0, a, ":")
        split(a[2], b, FS)
        gsub(",", "", b[1])
        tot = int(b[1])
    }
}

FNR != NR {
    sum += $1
}

END {
    print "Unique Reads: ", int(tot - sum)
    print "Duplicates:", sum
}
