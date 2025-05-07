idmap = open('name_to_array_id.txt','r')
id_dict = {}
for line in idmap:
    dogname,arrayid = line.strip().split(',')
    id_dict[arrayid] = dogname

vcfin = open('dogs_merged.vcf','r')
vcfout = open('namesfixed_dogs_merged.vcf','w')
for line in vcfin:
    if 'FORMAT' not in line:
        vcfout.write(line)
    else:
        linelist = line.strip().split()
        header = '\t'.join(linelist[:9])
        array_ids = linelist[9:]
        name_list = []
        for id in array_ids:
            name_list.append(id_dict[id])
        vcfout.write('{}\t{}\n'.format(header,'\t'.join(name_list)))

vcfout.close()
