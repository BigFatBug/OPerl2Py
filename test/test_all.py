# coding=utf8

from frb_res_fil import FrbResHandler
from com_mer_str import ComMerStrHandler
from snp_split import SnpSplitHandler
from class_pos_ano import ClassPosAnoHandler

def test_all():

    ori_vcf = '161125_IonXpress_077.freebayes.vcf'
    sam = '161125_IonXpress_077.bwa.quafil.idf.eff.sam'
    bed = 'MP33_lib.bed'
    fa = 'MP33_seq.fa'
    xls = 'mp33_pos_rev_foruse.xls'

    frb_handler = FrbResHandler()
    com_handler = ComMerStrHandler(bed, fa, sam)
    snp_handler = SnpSplitHandler()
    ano_handler = ClassPosAnoHandler(fa, xls)

    # # frb
    # read title
    frb_handler.read_title(ori_vcf)

    # exectube
    frb_handler.exectube(frb_handler.read_file(ori_vcf))

    # write file
    frb_handler.write_res_file('temp/out.freebayes.fil.res')
    frb_handler.write_vcf_out_file('temp/out.freebayes.fil.vcf')

    # # com
    # read title
    com_handler.read_title('temp/out.freebayes.fil.vcf')

    # com1_2
    com_handler.com1(com_handler.read_vcf('temp/out.freebayes.fil.vcf'))
    # com2_1
    com_handler.com2(com_handler.out1)
    # com3_2
    com_handler.com3(com_handler.out2)
    # com replace
    com_handler.com4(com_handler.out3, frb_handler.res_out)
    # get vcf out
    com_handler.filter_vcf(frb_handler.vcf_out_result)

    # write file
    with open('temp/out.sus.res', 'w') as a:
        for data in com_handler.out1:
            a.write(data + '\n')
    with open('temp/out.merge.res', 'w') as a:
        for data in com_handler.out2:
            a.write(data + '\n')
    with open('temp/out.check.res', 'w') as a:
        for data in com_handler.out3:
            a.write(data + '\n')
    with open('temp/out.mergerev.xls', 'w') as a:
        for data in com_handler.out4:
            a.write(data + '\n')
    with open('temp/out.mergerev.vcf', 'w') as a:
        for title in com_handler.vcf_title:
            a.write(title)
        for data in com_handler.vcf_out:
            a.write(data + '\n')

    # # snp (need com)
    # read title

    # snp
    snp_handler.snp1(com_handler.out4)
    # snp split
    snp_handler.snp2(snp_handler.out1)
    # get vcf out
    snp_handler.transfer_vcf(snp_handler.out2)

    # write file
    with open('temp/out.mergerev.snpsplit.xls', 'w') as a:
        for data in snp_handler.out1:
            a.write(data + '\n')
    with open('temp/out.mergerev.snpsplitrev.xls', 'w') as a:
        for data in snp_handler.out2:
            a.write(data + '\n')
    with open('temp/out.mergerev.snpsplitrev.vcf', 'w') as a:
        for title in com_handler.vcf_title:
            a.write(title)
        for data in snp_handler.vcf_out:
            a.write(data + '\n')

    # # ano
    # read title
    ano_handler.read_title('temp/out.mergerev.snpsplitrev.vcf')
    # class pos ano
    ano_handler.ano1(ano_handler.read_vcf('temp/out.mergerev.snpsplitrev.vcf'))
    # get vcf out
    ano_handler.transfer_vcf(ano_handler.out)
    # write file
    with open('temp/out.ano.xls', 'w') as a:
        for data in ano_handler.out:
            a.write(data + '\n')
    with open('temp/out.ano.vcf', 'w') as a:
        for title in ano_handler.vcf_title:
            a.write(title)
        for data in ano_handler.vcf_out:
            a.write(data + '\n')

    assert frb_handler.res_out
    assert frb_handler.vcf_out_result
    assert com_handler.out1
    assert com_handler.out2
    assert com_handler.out3
    assert com_handler.out4
    assert com_handler.vcf_out
    assert snp_handler.out1
    assert snp_handler.out2
    assert snp_handler.vcf_out
    assert ano_handler.out
    assert ano_handler.vcf_out

def test1():
    py_out = open('temp/out.freebayes.fil.res', 'r')
    perl_out = open('perl_out/1.freebayes.fil.res', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()

def test2():
    py_out = open('temp/out.sus.res', 'r')
    perl_out = open('perl_out/2.sus.res', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()

def test3():
    py_out = open('temp/out.merge.res', 'r')
    perl_out = open('perl_out/3.merge.res', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()

def test4():
    py_out = open('temp/out.check.res', 'r')
    perl_out = open('perl_out/4.check.res', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()

def test5():
    py_out = open('temp/out.mergerev.xls', 'r')
    perl_out = open('perl_out/5.mergerev.xls', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()

def test6():
    py_out = open('temp/out.mergerev.snpsplit.xls', 'r')
    perl_out = open('perl_out/6.mergerev.snpsplit.xls', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()

def test7():
    py_out = open('temp/out.mergerev.snpsplitrev.xls', 'r')
    perl_out = open('perl_out/7.mergerev.snpsplitrev.xls', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()

def test8():
    py_out = open('temp/out.ano.xls', 'r')
    perl_out = open('perl_out/8.ano.xls', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()