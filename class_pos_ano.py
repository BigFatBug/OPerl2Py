# coding=utf8
import sys
import HTSeq
import itertools

class ClassPosAnoHandler:

    def __init__(self, faFile, xlsFile):
        self.out = []
        self.vcf_out = []
        self.faFile = faFile
        self.xlsFile = xlsFile
        self.vcf_title = []

    def read_title(self, inputFile):
        with open(inputFile, 'r') as file:
            for data in file:
                if data[0] != '#':
                    return
                self.vcf_title.append(data)

    def read_vcf(self, inputFile):
        return [data for data in itertools.islice(HTSeq.VCF_Reader(inputFile), sys.maxsize)]

    def ano1(self, datas):
        cot = {}
        spr = {}
        for data in datas:
            cot[data.chrom][data.pos.start] = cot.setdefault(data.chrom, {}).setdefault(data.pos.start, 0) + 1
            # HA = 0
            info_split = data.info.split(';')
            if len(info_split) < 7:
                spr.setdefault(data.chrom, {}).setdefault(data.pos.start, {})[cot[data.chrom][data.pos.start]] = '%s\t%s\t%s\t%s' % (data.ref, data.alt[0], info_split[0].split('=')[1], info_split[1].split('=')[1])
            else:
                spr.setdefault(data.chrom, {}).setdefault(data.pos.start, {})[
                    cot[data.chrom][data.pos.start]] = '%s\t%s\t%s\t%s\t%s-%s-%s-%s' % (
                data.ref, data.alt[0], info_split[0].split('=')[1], info_split[1].split('=')[1], info_split[2].split('=')[1], info_split[3].split('=')[1],info_split[4].split('=')[1], info_split[5].split('=')[1])

        rfb = {}
        for fa_data in itertools.islice(HTSeq.FastaReader(self.faFile), sys.maxsize):
            inf = fa_data.name.split(':')
            for i in range(0, len(fa_data.seq)):
                rfb.setdefault(inf[0], {})[int(inf[1]) + i] = chr(fa_data.seq[i])

        exi = {}
        with open(self.xlsFile, 'r') as xls_reader:
            for xls_data in xls_reader.readlines():
                spl = xls_data[:-1].split('\t')

                if len(spl[3]) == len(spl[4]):
                    if cot.get(spl[0], {}).get(int(spl[1])):
                        for m in range(1, cot[spl[0]][int(spl[1])] + 1):
                            mid = spr[spl[0]][int(spl[1])][m].split('\t')
                            if mid[0] == spl[3] and mid[1] == spl[4]:
                                out = 'Class_1\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (spl[0], spl[5], spl[6], spl[7], spl[8], mid[2], mid[3])
                                if len(mid) > 4:
                                    out += '\t%s' % mid[4]
                                self.out.append(out)
                                exi.setdefault(spl[0], {}).setdefault(int(spl[1]), {}).setdefault(mid[0], {})[mid[1]] = 1
                                break

                else:
                    for i in range(int(spl[1]), int(spl[2]) + 1):
                        if cot.get(spl[0], {}).get(i):
                            for m in range(1, cot[spl[0]][i] + 1):
                                mid = spr[spl[0]][i][m].split('\t')
                                if len(mid[0]) == len(mid[1]) or len(mid[0]) + i - 1 > int(spl[2]):
                                    continue
                                lsq = ''
                                rsq = ''
                                msp = i + len(mid[0]) - 1

                                for n in range(int(spl[1]), i):
                                    lsq += rfb[spl[0]][n]
                                for n in range(msp + 1, int(spl[2]) + 1):
                                    rsq += rfb[spl[0]][n]

                                nfq = lsq + mid[0] + rsq
                                naq = lsq + mid[1] + rsq
                                if nfq == spl[3] and naq == spl[4]:
                                    out = 'Class_1\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (spl[0], spl[5], spl[6], spl[7], spl[8], mid[2], mid[3])
                                    if len(mid) > 4:
                                        out += '\t%s' % mid[4]
                                    self.out.append(out)
                                    exi.setdefault(spl[0], {}).setdefault(i, {}).setdefault(mid[0], {})[mid[1]] = 1
                                    break

        for data in datas:
            if not exi.get(data.chrom, {}).get(data.pos.start, {}).get(data.ref, {}).get(data.alt[0]):
                spp = data.pos.start + len(data.ref) - 1
                if 'HA=0' in data.info:
                    info_split = data.info.split(';')
                    out = 'Undefined\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (data.chrom, data.pos.start, spp, data.ref, data.alt[0], info_split[0].split('=')[1], info_split[1].split('=')[1])
                    out += '\t%s-%s-%s-%s' % (info_split[2].split('=')[1], info_split[3].split('=')[1],info_split[4].split('=')[1], info_split[5].split('=')[1])
                    self.out.append(out)

    def transfer_vcf(self, datas):
        for data in datas:
            d = data.split('\t')
            info8 = d[8].split('-')
            info = 'DP=%s;AO=%s;FS=%s;RS=%s;SAF=%s;SAR=%s;HA=0;Level=%s' % (d[6], d[7], info8[0],info8[1],info8[2],info8[3], d[0])
            self.vcf_out.append('\t'.join([d[1], d[2], '.', d[4], d[5], '.', '.', info]))

if __name__ == '__main__':
    c = ClassPosAnoHandler('test/MP33_seq.fa', 'test/mp33_pos_rev_foruse.xls')
    c.ano1(c.read_vcf('test/temp/file9'))
    c.transfer_vcf()
    with open('over', 'w') as a:
        for x in c.out:
            a.write(x + '\n')