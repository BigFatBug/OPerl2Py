# coding=utf8
import sys
import copy, re
import HTSeq
import itertools

class FrbResHandler:

    def __init__(self, dp_replace_flag=True):
        self.vcf_title = []
        self.vcf_out_result = []
        self.res_out = []
        self.filter_index = []
        self.temp_vcf_data = None
        self.dp_replace_flag = dp_replace_flag

    def read_title(self, inputFile):
        with open(inputFile, 'r') as file:
            for data in file:
                if data[0] != '#':
                    return
                self.vcf_title.append(data)

    def read_file(self, inputFile):
        return [data for data in itertools.islice(HTSeq.VCF_Reader(inputFile), sys.maxsize)]

    def exectube(self, datas):
        for data in datas:
            ift, _, inf = data.info.split('\t')
            inf = inf.split(':')
            ift = ift.split(';')
            bdp = [int(i) for i in inf[2].split(',')]
            qus = [int(i) for i in inf[6].split(',')]
            msum = 0
            for i in bdp:
                msum += int(i)
            self.temp_vcf_data = self.split_vcf_data(data, msum)
            afs = '0'
            ars = '0'
            rff = 0
            rrr = 0
            con = 0
            ath = {}

            for i in ift:
                if 'SAF' in i:
                    afs = i.replace('SAF=', '')
                if 'SAR' in i:
                    ars = i.replace('SAR=', '')
                if 'SRF' in i:
                    rff = int(i.replace('SRF=', ''))
                if 'SRR' in i:
                    rrr = int(i.replace('SRR=', ''))

            aff = afs.split(',')
            for i in aff:
                rff += int(i)

            arr = ars.split(',')
            for i in arr:
                rrr += int(i)

            if len(data.ref) < 3:
                self.write(data, msum, bdp, qus, rff, rrr, aff, arr)

            elif len(data.ref) == 3:
                rfb = data.ref
                if rfb[2] != rfb[1] and rfb[1] != rfb[0]:
                    self.write(data, msum, bdp, qus, rff, rrr, aff, arr)
                else:
                    if rfb[2] == rfb[1]:
                        midr = rfb + rfb[2]
                        midd = rfb
                        if midd[-1] in ['A', 'T', 'C', 'G']:
                            midd = midd[:-1]
                    else:
                        midr = rfb + data.ref
                        midd = data.ref
                        if midd[0] in ['A', 'T', 'C', 'G']:
                            midd = midd[1:]
                    self.deal1(midd, data, ath, msum, bdp, midr, qus, rff, rrr, aff, arr, con)

            else:
                rfb = data.ref
                if rfb[-1] == rfb[-2]:
                    ttb = rfb[1]
                    btl = 1
                    for i in range(len(rfb)):
                        btl = 0 if rfb[i] != ttb else btl
                    if btl:
                        midr = rfb + ttb
                        midd = rfb
                        if midd[-1] in ['A', 'T', 'C', 'G']:
                            midd = midd[:-1]
                        self.deal1(midd, data, ath, msum, bdp, midr, qus, rff, rrr, aff, arr, con)
                    else:
                        self.write(data, msum, bdp, qus, rff, rrr, aff, arr)

                else:
                    ttb = rfb[1]
                    btl = 1
                    for i in range(1, len(rfb) - 1):
                        btl = 0 if rfb[i] != ttb else btl
                    if btl:
                        midr = rfb
                        midd = rfb
                        if midr[-1] in ['A', 'T', 'C', 'G']:
                            midr = midr[:-1]
                        if midd[-1] in ['A', 'T', 'C', 'G'] and midd[-2] in ['A', 'T', 'C', 'G']:
                            midd = midd[:-2]
                        midr = midr + ttb
                        self.deal2(midd, data, ath, msum, bdp, midr, qus, rff, rrr, aff, arr, con)

                    else:
                        self.write(data, msum, bdp, qus, rff, rrr, aff, arr)

    def write(self, data, msum, bdp, qus, rff, rrr, aff, arr):
        self.extend()
        if len(data.alt) > 1:
            att = data.alt
            for i in range(len(att)):
                self.write1(data, msum, bdp, qus, rff, rrr, aff, arr, att, i)
        else:
            self.write2(data, msum, bdp, qus, rff, rrr, aff, arr)

    def write1(self, data, msum, bdp, qus, rff, rrr, aff, arr, att, i):
        self.res_out.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
            data.chrom, data.pos.start, data.ref, att[i], msum, bdp[i+1], qus[i], rff ,rrr, aff[i], arr[i]))

    def write2(self, data, msum, bdp, qus, rff, rrr, aff, arr):
        self.res_out.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
            data.chrom, data.pos.start, data.ref, ','.join(data.alt), msum, bdp[1], qus[0], rff, rrr, aff[0], arr[0]))

    def deal1(self, midd, data, ath, msum, bdp, midr, qus, rff, rrr, aff, arr, con):
        if len(data.alt) > 1:
            att = data.alt
            for i in range(len(att)):
                ath[i+1] = att[i]
                con = i + 1
            for i in range(1, len(att) + 1):
                if ath[i] in midd or midr in ath[i]:
                    bdp[0] += bdp[i]
                    ath.pop(i)
                    con -= 1
            if con:
                for i in range(0, len(att)):
                    if ath.get(i + 1):
                        self.write1(data, msum, bdp, qus, rff, rrr, aff, arr, att, i)
                    else:
                        self.temp_vcf_data[i] = None
                self.extend()
        else:
            if data.alt[0] in midd or midr in data.alt[0]:
                if bdp[1] > bdp[0]:
                    self.extend()
                    self.write2(data, msum, bdp, qus, rff, rrr, aff, arr)
            else:
                self.extend()
                self.write2(data, msum, bdp, qus, rff, rrr, aff, arr)

    def deal2(self, midd, data, ath, msum, bdp, midr, qus, rff, rrr, aff, arr, con):
        if len(data.alt) > 1:
            att = data.alt
            for i in range(len(att)):
                ath[i + 1] = att[i]
                if ath[i+1][-1] in ['A', 'T', 'C', 'G']:
                    ath[i + 1] = ath[i+1][:-1]
                con = i + 1
            for i in range(1, len(att) + 1):
                if ath[i] in midd or midr in ath[i]:
                    bdp[0] += bdp[i]
                    ath.pop(i)
                    con -= 1
            if con:
                for i in range(0, len(att)):
                    if ath.get(i + 1):
                        self.write1(data, msum, bdp, qus, rff, rrr, aff, arr, att, i)
                    else:
                        self.temp_vcf_data[i] = None
                self.extend()
        else:
            athm = data.alt[0]
            if athm[-1] in ['A', 'T', 'C', 'G']:
                athm = athm[:-1]
            if athm in midd or midr in athm:
                if bdp[1] > bdp[0]:
                    self.extend()
                    self.write2(data, msum, bdp, qus, rff, rrr, aff, arr)
            else:
                self.extend()
                self.write2(data, msum, bdp, qus, rff, rrr, aff, arr)

    def split_vcf_data(self, data, msum):
        info = str(data._original_line[:-1])
        if self.dp_replace_flag:
            info = re.sub(';DP=\d+;', ';DP=%s;' % msum, info)
        result = []
        if len(data.alt) <= 1:
            result.append(info)
        else:
            details = info.split('\t')
            line8_split = [r.split(',') for r in details[7].split(';')]
            line10_datas = []
            line10_split = details[9].split(':')
            for i in range(len(data.alt)):
                line10_data = copy.deepcopy(line10_split)
                for j in range(len(line10_data)):
                    x = line10_data[j].split(',')
                    if len(x) == len(data.alt) + 1:
                        line10_data[j] = x[0] + ',' + x[i+1]
                    elif len(x) == len(data.alt):
                        line10_data[j] = x[i]
                    elif len(x) > len(data.alt):
                        line10_data[j] = ','.join([x[k] for k in range(len(x)) if (k+i) % len(data.alt) == 0])
                line10_datas.append(':'.join(line10_data))

            for i in range(len(data.alt)):
                row_data = copy.deepcopy(details)
                row_data[4] = data.alt[i]
                row_data[7] = ';'.join([((re.sub('=[0-9a-zA-Z\.]+', '='+l[i], l[0]) if i else l[i]) if len(l) > 1 else l[0]) for l in line8_split])

                row_data[9] = line10_datas[i]
                result.append('\t'.join(row_data))
        return result

    def write_vcf_out_file(self, vcfOutFile):
        with open(vcfOutFile, 'w') as out:
            self.write_title(out)
            for i in self.vcf_out_result:
                out.write(i + '\n')

    def write_res_file(self, resFile):
        with open(resFile, 'w') as out:
            for i in self.res_out:
                out.write(i + '\n')

    def write_title(self, out):
        for data in self.vcf_title:
            out.write(data)

    def extend(self):
        for t in self.temp_vcf_data:
            if t:
                self.vcf_out_result.append(t)

if __name__ == '__main__':
    frb = FrbResHandler()
    frb.read_title('test/160822_IonXpress_3.freebayes.vcf')
    frb.exectube(frb.read_file('test/160822_IonXpress_3.freebayes.vcf'))
    frb.write_vcf_out_file('1')
    frb.write_res_file('2')