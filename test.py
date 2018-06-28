import json

count_total_reads = 0
count_unmapped = 0
count_reg = 0
count_norm = 0
count_collisions = 0
count_lowqcollisions = 0
count_mapq0 = 0
for norm in norm_res:
    with open(norm,'r') as fp:
        values = fp.read().split()
        count_total_reads += values[0]
        count_unmapped += values[1]
        count_reg += values[2]
        count_norm += values[3]
        count_collisions += values[4]
        count_lowqcollisions += values[5]
        count_mapq0 += values[6]
    fp.close()
data = {"Total reads": count_total_reads, "Total Unmapped": count_unmapped, "Total regular": count_reg,"Total normal": count_norm, "Total collisions": count_collisions, "Total lowqcollisions": count_lowqcollisions,"Total mapq0": count_mapq0}

with open('align_qc.json', 'w') as outfile:
    json.dump(data, outfile)
outfile.close()


