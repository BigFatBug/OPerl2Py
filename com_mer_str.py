# coding=utf8
import sys
import re
import HTSeq
import itertools

class ComMerStrHandler:

    """
    inputFile: *.res
    configFile: *.bed
    faFile: *.fa
    samFile: *.sam
    """
    def __init__(self, bedFile, faFile, samFile):
        self.bedFile = bedFile
        self.faFile = faFile
        self.samFile = samFile

        self.out1 = []
        self.out2 = []
        self.out3 = []
        self.out4 = []

        self.vcf_out = []
        self.vcf_title = []

    def read_title(self, inputFile):
        with open(inputFile, 'r') as file:
            for data in file:
                if data[0] != '#':
                    return
                self.vcf_title.append(data)

    def read_vcf(self, inputFile):
        return [data for data in itertools.islice(HTSeq.VCF_Reader(inputFile), sys.maxsize)]

    def com1(self, datas):
        con = {}
        frt = {}
        frn = {}
        cnt = {}
        tct = []
        for data in datas:
            if 'POLY=1' in data._original_line:
                continue
            dp = re.match('.*;DP=(\d+)', data.info).groups()[0]
            ao = re.match('.*;AO=(\d+)', data.info).groups()[0]
            t = [data.chrom, data.pos.start, data.ref, data.alt[0], dp, ao]
            con[t[0]][t[1]] = con.setdefault(t[0], {}).setdefault(t[1], 0) + 1
            frt.setdefault(t[0], {}).setdefault(t[1], {})[con[t[0]][t[1]]] = '%s-%s' % (t[2], t[3])
            frn.setdefault(t[0], {}).setdefault(t[1], {})[con[t[0]][t[1]]] = '%s-%s' % (t[4], t[5])

        for deb_data in itertools.islice(HTSeq.BED_Reader(self.bedFile), sys.maxsize):
            for i in range(deb_data.iv.start, deb_data.iv.end):
                if con.get(deb_data.iv.chrom, {}).get(i):
                    bod = min(i+50, deb_data.iv.end)
                    for m in range(1, con[deb_data.iv.chrom][i] + 1):
                        res = '%s-%s-%s' % (deb_data.iv.chrom, i, frt[deb_data.iv.chrom][i][m])
                        lef = frn[deb_data.iv.chrom][i][m].split('-')
                        mko = '%s-%s' % (i, frt[deb_data.iv.chrom][i][m])
                        rfn = int(lef[0]) - int(lef[1])
                        atn = int(lef[1])
                        tel = 0
                        for n in range(i+1, bod+1):
                            if con.get(deb_data.iv.chrom, {}).get(n):
                                for a in range(1, con[deb_data.iv.chrom][n]+1):
                                    mkt = '%s-%s' % (n, frt[deb_data.iv.chrom][n][a])
                                    if not cnt.get(mko, {}).get(mkt):
                                        mid = frn[deb_data.iv.chrom][n][a].split('-')
                                        rft = int(mid[0]) - int(mid[1])
                                        att = int(mid[1])
                                        kfz = (rfn + atn + rft + att) * (rfn * att - rft * atn) ** 2 / ((rfn + atn) * (rft + att) * (rfn + rft) * (atn + att));
                                        if kfz < 2.5:
                                            res = '%s\t%s-%s-%s' % (res, deb_data.iv.chrom, n, frt[deb_data.iv.chrom][n][a])
                                            tct.append(mkt)
                                            rfn += rft
                                            atn += att
                                            tel = 1
                                            break
                        if tel:
                            for n in range(0, len(tct)):
                                for a in range(n+1, len(tct)):
                                    cnt.setdefault(tct[n], {})[tct[a]] = 1
                            self.out1.append('%s\t%s\t%s' % (rfn, atn, res))

    def com2(self, datas):
        mrq = {}
        maq = {}
        rfq = {}
        for fa_data in itertools.islice(HTSeq.FastaReader(self.faFile), sys.maxsize):
            inf = fa_data.name.split(':')
            for i in range(0, len(fa_data.seq)):
                rfq.setdefault(inf[0], {})[int(inf[1]) + i] = chr(fa_data.seq[i])

        for data in datas:
            out = data.replace('\n', '')
            spl = out.split('\t')
            seqr = ''
            seqa = ''
            cho = ''
            lef = 100000000000000000
            rig = 0
            for i in range(2, len(spl)):
                mid = spl[i].split('-')
                cho = mid[0]
                lef = min(lef, int(mid[1]))
                rig = max(rig, int(mid[1]))
                mrq.setdefault(cho, {})[int(mid[1])] = mid[2]
                maq.setdefault(cho, {})[int(mid[1])] = mid[3]

            it = lef
            while it <= rig:
                if mrq.get(cho, {}).get(it):
                    seqr += mrq[cho][it]
                    seqa += maq[cho][it]
                    it += len(mrq[cho][it])
                else:
                    seqr += str(rfq[cho][it])
                    seqa += str(rfq[cho][it])
                    it += 1

            rig = lef + len(seqr) - 1
            if len(seqr) >= 40:
                while seqr and seqa and seqr[0] == seqa[0]:
                    seqr = seqr[1:]
                    seqa = seqa[1:]
                    lef += 1
                while seqr and seqa and seqr[-1] == seqa[-1]:
                    seqr = seqr[:-1]
                    seqa = seqa[:-1]
                    rig -= 1
            self.out2.append('%s\t%s\t%s\t%s\t%s\t%s' % (cho, lef, rig, seqr, seqa, out))

    def com3(self, datas):
        sam_datas = [sam_data for sam_data in itertools.islice(HTSeq.SAM_Reader(self.samFile), sys.maxsize)]
        for data in datas:
            out = data.replace('\n', '')
            eif = out.split('\t')
            exchr = eif[0]
            exlef = int(eif[1])
            exrig = int(eif[2])
            rfseq = eif[3]
            atseq = eif[4]
            exrfn = int(eif[5])
            exatn = int(eif[6])

            atn = 0
            rfn = 0
            rfn_ops = 0
            rfn_nes = 0
            atn_ops = 0
            atn_nes = 0

            for sam_data in sam_datas:
                ppt = sam_data.iv.start
                if sam_data.iv.chrom == exchr:
                    seq = ''
                    cig = ''
                    for cigar in sam_data.cigar:
                        if cigar.type in 'MIDSH':
                            clen = cigar.size
                            typ = cigar.type
                            if typ != 'H':
                                cig += clen * typ
                    bas = str(sam_data.read_as_aligned)
                    m = 0
                    you = 0
                    sto = 1 if ppt < exlef else 0
                    for i in range(0, len(cig)):
                        if cig[i] == 'S':
                            m += 1
                        elif cig[i] == 'M':
                            ppt += 1
                            if ppt >= exlef and ppt <= exrig:
                                seq += bas[m]
                                you = 1
                            m += 1
                        elif cig[i] == 'I':
                            if ppt >= exlef and ppt <= exrig:
                                seq += bas[m]
                                you = 1
                            m += 1
                        elif cig[i] == 'D':
                            ppt += 1
                    sto += 1 if ppt >= exrig else 0

                    if not you:
                        continue

                    if sto == 2:
                        if seq == atseq:
                            atn += 1
                            if sam_data.flag == 0:
                                atn_ops += 1
                            else:
                                atn_nes += 1
                        else:
                            if len(rfseq) == len(atseq):
                                amn, rmn = 0, 0
                                for n in range(0, len(rfseq)):
                                    if n >= len(seq):
                                        amn += 1
                                        rmn += 1
                                    else:
                                        amn += 1 if atseq[n] != seq[n] else 0
                                        rmn += 1 if rfseq[n] != seq[n] else 0
                                if amn < rmn:
                                    atn += 1
                                    if sam_data.flag == 0:
                                        atn_ops += 1
                                    else:
                                        atn_nes += 1
                                else:
                                    rfn += 1
                                    if sam_data.flag == 0:
                                        rfn_ops += 1
                                    else:
                                        rfn_nes += 1
                            else:
                                if abs(len(atseq) - len(seq)) < abs(len(rfseq) - len(seq)):
                                    atn += 1
                                    if sam_data.flag == 0:
                                        atn_ops += 1
                                    else:
                                        atn_nes += 1
                                else:
                                    rfn += 1
                                    if sam_data.flag == 0:
                                        rfn_ops += 1
                                    else:
                                        rfn_nes += 1

            llb = (rfn + atn + exrfn + exatn) * (rfn * exatn - atn * exrfn) ** 2 / ((rfn + atn) * (exrfn + exatn) * (rfn + exrfn) * (atn + exatn))
            if llb < 3.841:
                ops = rfn_ops + atn_ops
                nes = rfn_nes + atn_nes
                self.out3.append('%s-%s\t%s-%s-%s-%s\t%s' % (rfn, atn, ops, nes, atn_ops, atn_nes, out))

    def com4(self, datas, ori_datas):
        alt = {}
        exi = {}
        for data in datas:
            spl = data.split('\t')
            ton = spl[0].split('-')
            dep = int(ton[0]) + int(ton[1])
            spl[1] = spl[1].replace('-', '\t')
            mid = '%s\t%s\t%s\t%s\t%s\t%s\t512\t%s' % (spl[2], spl[3], spl[5], spl[6], dep, ton[1], spl[1])
            for i in range(9, len(spl)):
                tbk = spl[i].split('-')
                alt.setdefault(tbk[0], {}).setdefault(tbk[1], {}).setdefault(tbk[2], {})[tbk[3]] = mid

        for ori_data in ori_datas:
            ori_data = ori_data.replace('\n', '')
            spl = ori_data.split('\t')
            if alt.get(spl[0], {}).get(spl[1], {}).get(spl[2], {}).get(spl[3]):
                if exi.get(alt[spl[0]][spl[1]][spl[2]][spl[3]]):
                    self.out4.append('%s\thahaha' % ori_data)
                    continue
                else:
                    self.out4.append(alt[spl[0]][spl[1]][spl[2]][spl[3]] + '\t%s' % ori_data.split('\t')[-1])
                    exi[alt[spl[0]][spl[1]][spl[2]][spl[3]]] = 1
                    self.out4.append('%s\thahaha' % ori_data)
                    continue
            else:
                if 15 * int(spl[5]) <= int(spl[6]):
                    self.out4.append(ori_data)

    def filter_vcf(self, frb_vcf_out):
        i = 0
        j = 0
        while i < len(self.out4) and j < len(frb_vcf_out):
            i_split = self.out4[i].split('\t')
            j_split = frb_vcf_out[j].split('\t')
            if i_split[0] + i_split[1] == j_split[0] + j_split[1]:
                self.vcf_out.append(frb_vcf_out[j])
                j += 1
            else:
                i += 1