# Function to split Vcf record and return Ordered dict
def vcfRecord(line, info_format):

    splittedLine = line.split('\t')

    entry = {
        "CHROM": splittedLine[0],
        "POS": splittedLine[1],
        "ID": splittedLine[2],
        "REF": splittedLine[3],
        "ALT": splittedLine[4],
        "QUAL": splittedLine[5],
        "FILTER": splittedLine[6].split(';'),
        "INFO": {
            i.split('=')[0]: i.split('=')[1].split(',')
            for i in splittedLine[7].split(';')
        },
        "FORMAT": {
            "GT": splittedLine[9].split(':')[0],
            "VAF": splittedLine[9].split(':')[1],
            "VD": splittedLine[9].split(':')[2],
            "DP": splittedLine[9].split(':')[3]
        }
    }

    info_csq_reformatted = {}
    for i in info_format:
        info_csq_reformatted[i] = []
        for j in entry['INFO']['CSQ']:
            info_csq_reformatted[i].append(j.split('|')[info_format.index(i)])

    entry['INFO']['CSQ'] = info_csq_reformatted

    return entry, len(info_csq_reformatted[info_format[0]])

#Fuction to get rsids from the ordered dict
def get_rsIDs(existing_variation):
    rsids = []
    for r in existing_variation:
        if r != '':
            rs = r.split('&')
            for rs_element in rs:
                if rs_element not in rsids and rs_element.startswith('rs'):
                    rsids.append(rs_element)
    return rsids

#Fucntion to add rsIds to the ID column of the VCF
def reconstructAnnotation_withrsID(vcfEntry, transcriptCounts,rsid_list=[]):

    variant_callers = f"variant_callers={'|'.join(vcfEntry['INFO']['variant_callers'])}"

    reconstruted_CSQList = []
    z = 0
    while z < transcriptCounts:
        temp_CSQ_element = '|'.join(vcfEntry['INFO']['CSQ'][k][z]
                                for k in vcfEntry['INFO']['CSQ'].keys())
        reconstruted_CSQList.append(temp_CSQ_element)
        z += 1

    reconstruted_CSQRecord = f"CSQ={','.join(reconstruted_CSQList)}"
    reconstruted_INFORecord = f"{variant_callers};{reconstruted_CSQRecord}"

    newVCFEntryList = []

    if rsid_list == []:
        newVCFEntryList.append(f"{vcfEntry['CHROM']}\t{vcfEntry['POS']}\t{vcfEntry['ID']}\t{vcfEntry['REF']}\t{vcfEntry['ALT']}\t{vcfEntry['QUAL']}\t{'|'.join(vcfEntry['FILTER'])}\t{reconstruted_INFORecord}\tGT:VAF:VD:DP\t{':'.join(map(str, vcfEntry['FORMAT'].values()))}")

    else:
        for rsID in rsid_list:
            newVCFEntryList.append(
                f"{vcfEntry['CHROM']}\t{vcfEntry['POS']}\t{rsID}\t{vcfEntry['REF']}\t{vcfEntry['ALT']}\t{vcfEntry['QUAL']}\t{'|'.join(vcfEntry['FILTER'])}\t{reconstruted_INFORecord}\tGT:VAF:VD:DP\t{':'.join(map(str, vcfEntry['FORMAT'].values()))}"
            )

    newVCFEntry = '\n'.join(newVCFEntryList)

    return newVCFEntry
