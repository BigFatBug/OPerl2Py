# coding=utf8
import sys
import HTSeq
import itertools

class SnpSplitHandler:

    def __init__(self):
        self.out1 = []
        self.out2 = []
        self.vcf_out = []

    def get_out_data(self, spl, index):
        if index == 1:
            return '%s\t%s\t%s\t%s\t%s\t%s\thahaha' % (spl[0], spl[1], spl[2], spl[3], spl[4], spl[5])
        if index == 2:
            return '%s\t%s\t%s\t%s\t%s\t%s' % (spl[0], spl[1], spl[2], spl[3], spl[4], spl[5])
        if index == 3:
            return '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\thahaha' % (spl[0], spl[1], spl[2], spl[3], spl[4], spl[5], spl[7], spl[8], spl[9], spl[10])
        if index == 4:
            return '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (spl[0], spl[1], spl[2], spl[3], spl[4], spl[5], spl[7], spl[8], spl[9], spl[10])

    def snp1(self, datas):
        for data in datas:
            spl = data.split('\t')
            if len(spl[2]) != len(spl[3]):
                if len(spl) <= 10:
                    if spl[-1] == 'hahaha':
                        self.out1.append(self.get_out_data(spl, 1))
                    else:
                        self.out1.append(self.get_out_data(spl, 2))
                else:
                    if spl[-1] == 'hahaha':
                        self.out1.append(self.get_out_data(spl, 3))
                    else:
                        self.out1.append(self.get_out_data(spl, 4))
            else:
                san = 0
                for i in range(0, len(spl[2])):
                    san += 1 if spl[2][i] == spl[3][i] else 0
                if 2 * san >= len(spl[2]):
                    tot = 0
                    pot = {}
                    bsr = {}
                    bsa = {}
                    rba = {}
                    cnt = {}

                    for i in range(0, len(spl[2])):
                        rba[int(spl[1]) + i] = spl[2][i]
                    for i in range(0, len(spl[2])):
                        if spl[2][i] != spl[3][i]:
                            tot += 1
                            pot[tot] = int(spl[1]) + i
                            bsr[tot] = spl[2][i]
                            bsa[tot] = spl[3][i]
                    for i in range(1, tot):
                        if pot[i+1] - pot[i] < 3:
                            cnt[i] = i + 1
                    for i in range(1, tot):
                        if cnt.get(i):
                            ano = cnt[i]
                            while cnt.get(ano):
                                cnt[i] = cnt[ano]
                                mid = ano
                                ano = cnt[ano]
                                cnt.pop(mid)
                    for i in range(1, tot):
                        if cnt.get(i):
                            rsq = bsr[i]
                            asq = bsa[i]
                            for m in range(i+1, cnt[i]+1):
                                for n in range(pot.get(m-1, 0) + 1, pot.get(m, 0)):
                                    rsq += rba.get(n, '')
                                    asq += rba.get(n, '')
                                rsq += bsr[m]
                                asq += bsa[m]
                                pot.pop(m)
                            bsr[i] = rsq
                            bsa[i] = asq
                    for i in range(1, tot+1):
                        if pot.get(i):
                            if len(spl) <= 10:
                                self.out1.append('%s\t%s\t%s\t%s\t%s\t%s' % (spl[0], pot[i], bsr[i], bsa[i], spl[4], spl[5]))
                            else:
                                self.out1.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (spl[0], pot[i], bsr[i], bsa[i], spl[4], spl[5], spl[7], spl[8], spl[9], spl[10]))

                else:
                    if len(spl) <= 10:
                        if spl[-1] == 'hahaha':
                            self.out1.append(self.get_out_data(spl, 1))
                        else:
                            self.out1.append(self.get_out_data(spl, 2))
                    else:
                        if spl[-1] == 'hahaha':
                            self.out1.append(self.get_out_data(spl, 3))
                        else:
                            self.out1.append(self.get_out_data(spl, 4))

    def snp2(self, datas):
        con = {}
        tel = {}
        cal = {}
        tim = {}

        for data in datas:
            spl = data.split('\t')
            con[spl[0]][spl[1]][spl[2]][spl[3]] = con.setdefault(spl[0], {}).setdefault(spl[1], {}).setdefault(spl[2], {}).setdefault(spl[3], 0) + 1
            if spl[-1] == 'hahaha':
                tel.setdefault(spl[0], {}).setdefault(spl[1], {}).setdefault(spl[2], {})[spl[3]] = 1

        for data in datas:
            spl = data.split('\t')
            if con[spl[0]][spl[1]][spl[2]][spl[3]] > 1:
                if tel.get(spl[0], {}).get(spl[1], {}).get(spl[2], {}).get(spl[3]):
                    if spl[-1] == 'hahaha':
                        self.out2.append(data.replace('\thahaha', ''))
                    else:
                        continue
                else:
                    if tim.get(spl[0], {}).get(spl[1], {}).get(spl[2], {}).get(spl[3]):
                        mid = cal[spl[0]][spl[1]][spl[2]][spl[3]].split('-')
                        spl[5] = str(int(spl[5]) + int(mid[0]))
                        spl[8] = str(int(spl[8]) + int(mid[1]))
                        spl[9] = str(int(spl[9]) + int(mid[2]))
                        self.out2.append('\t'.join(spl))
                    else:
                        cal.setdefault(spl[0], {}).setdefault(spl[1], {}).setdefault(spl[2], {})[spl[3]] = '%s-%s-%s' % (spl[5], spl[8], spl[9])
                        tim.setdefault(spl[0], {}).setdefault(spl[1], {}).setdefault(spl[2], {})[spl[3]] = 1
            else:
                self.out2.append(data)

    def transfer_vcf(self, datas):
        for data in datas:
            d = data.split('\t')
            info = 'DP=%s;AO=%s;FS=%s;RS=%s;SAF=%s;SAR=%s;HA=%s' % (
            d[4], d[5], d[6], d[7], d[8], d[9], 1 if d[-1] == 'hahaha' else 0)
            self.vcf_out.append('\t'.join([d[0], d[1], '.', d[2], d[3], '.', '.', info]))