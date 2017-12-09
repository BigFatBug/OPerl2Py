# coding=utf8

from frb_res_fil import FrbResHandler

def test_run1():
    # 运行一次测试用例，从输出文件中抽样检查
    handler = FrbResHandler()
    handler.exectube(handler.read_file('160822_IonXpress_20.freebayes.vcf'))
    handler.write_vcf_out_file('test1.out.vcf')
    handler.write_res_file('test1.res')
    with open('test1.res', 'r') as readFile:
        datas = readFile.readlines()
        assert len(datas) == 72
        assert datas[3].split('\t') == ['chr10','43613843','G','T','2084','2084','71289','1475','609','1475','609\n']
        assert datas[3].split('\t') == ['chr10','43613843','G','T','2084','2084','71289','1475','609','1475','609\n']
        assert datas[9].split('\t')[0] == 'chr11'
        assert sum([int(data.split('\t')[0].replace('chr', '')) for data in datas]) == 742

def test_run2():
    # 运行一次测试用例，从输出文件中抽样检查
    handler = FrbResHandler()
    handler.exectube(handler.read_file('160923_IonXpress_003.freebayes.vcf'))
    handler.write_vcf_out_file('test2.out.vcf')
    handler.write_res_file('test2.res')
    with open('test2.res', 'r') as readFile:
        datas = readFile.readlines()
        assert len(datas) == 45
        assert datas[0].split('\t') == ['chr1','115256501','TTGGTCT','CTGGTCC','1344','24','839','712','632','16','8\n']
        assert datas[9].split('\t') == ['chr11','533533','TGC','TGGC','915','61','1752','338','577','61','0\n']
        assert datas[-1].split('\t')[0] == 'chr9'
        assert sum([int(data.split('\t')[0].replace('chr', '')) for data in datas]) == 485

if __name__ == '__main__':
    test_run2()
