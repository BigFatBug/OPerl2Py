# coding=utf8

from frb_res_fil import FrbResHandler
from com_mer_str import ComMerStrHandler
from snp_split import SnpSplitHandler

def test_all():

    vcf = '161125_IonXpress_077.freebayes.vcf'
    sam = '161125_IonXpress_077.bwa.quafil.idf.eff.sam'
    bed = 'MP33_lib.bed'
    fa = 'MP33_seq.fa'

    frb_handler = FrbResHandler()
    com_handler = ComMerStrHandler(bed, fa, sam)
    snp_handler = SnpSplitHandler()

    datas = frb_handler.read_file(vcf)
    frb_handler.read_title(vcf)
    frb_handler.exectube(datas)

    frb_handler.write_res_file('temp/file1')

    frb_handler.write_vcf_out_file('temp/file2')

    # com1_2
    com_handler.com1(com_handler.read_vcf('temp/file2'))
    # com2_1
    com_handler.com2(com_handler.out1)
    # # com3_2
    com_handler.com3(com_handler.out2)
    # # com replace
    com_handler.com4(com_handler.out3, frb_handler.res_out)
    # # snp
    snp_handler.snp1(com_handler.out4)
    # # snp split
    snp_handler.snp2(snp_handler.out1)

    with open('temp/3', 'w') as a:
        for data in com_handler.out1:
            a.write(data + '\n')
    with open('temp/file4', 'w') as a:
        for data in com_handler.out2:
            a.write(data + '\n')
    with open('temp/file5', 'w') as a:
        for data in com_handler.out3:
            a.write(data + '\n')
    with open('temp/file6', 'w') as a:
        for data in com_handler.out4:
            a.write(data + '\n')
    with open('temp/file7', 'w') as a:
        for data in snp_handler.out1:
            a.write(data + '\n')
    with open('temp/file8', 'w') as a:
        for data in snp_handler.out2:
            a.write(data + '\n')

    # com replace filter vcf
    com_handler.filter_vcf(frb_handler.vcf_out_result)
    with open('temp/com_replace_out.vcf', 'w') as a:
        # vcf title
        for title in frb_handler.vcf_title:
            a.write(title)

        for data in com_handler.vcf_out:
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

def test1():
    py_out = open('temp/file1', 'r')
    perl_out = open('perl_out/1.freebayes.fil.res', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()

def test2():
    py_out = open('temp/file3', 'r')
    perl_out = open('perl_out/2.sus.res', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()

def test3():
    py_out = open('temp/file4', 'r')
    perl_out = open('perl_out/3.merge.res', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()

def test4():
    py_out = open('temp/file5', 'r')
    perl_out = open('perl_out/4.check.res', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()

def test5():
    py_out = open('temp/file6', 'r')
    perl_out = open('perl_out/5.mergerev.xls', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()

def test6():
    py_out = open('temp/file7', 'r')
    perl_out = open('perl_out/6.mergerev.snpsplit.xls', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()

def test7():
    py_out = open('temp/file8', 'r')
    perl_out = open('perl_out/7.mergerev.snpsplitrev.xls', 'r')
    assert perl_out.readlines() == py_out.readlines()
    py_out.close()
    perl_out.close()