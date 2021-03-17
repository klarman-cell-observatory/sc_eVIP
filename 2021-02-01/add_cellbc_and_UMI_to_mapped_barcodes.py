import argparse
import gzip
from itertools import islice
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--mapped')
    parser.add_argument('--R1')
    parser.add_argument('--out')
    parser.add_argument('--cell_bc_start',type=int,default=0)
    parser.add_argument('--cell_bc_end',type=int,default=16)
    parser.add_argument('--UMI_start',type=int,default=16)
    parser.add_argument('--UMI_end',type=int,default=26)
    args = parser.parse_args()

    out=gzip.open(args.out,'w')
    out_cells=gzip.open(args.out+'.cells.gz','w')
    mapped_reads={}
    readlist=set()
    for line in gzip.open(args.mapped,'r'):
        items=line.strip().split('\t')
        #['E00173:665:HH2H3CCXY:1:1101:26697:39651 2:N:0:0', '+', 'TP53_H193L', '52', 'CTGACTCAGACTGAGACAAAATGCCCGGGCTGCAAGTTCGGAAGAAGCTAGCGTGGTGTTAATTAACCTAGAGGGCCCGTTTAAACCCGCTGATCAGC', '-AAFFJFJJJAJ<----<F--FF--AFAA7AF7<JJFJ-FAFA-777AFJJFFJ--7F--F<FJJAJ77<J-AAFAJJJJF7J--F7<JA<J7<AAA-', '0']
        readname,strand,mapped_name,pos0,read,quality,num_alignments=items[0],items[1],items[2],items[3],items[4],items[5],items[6]
        readname=readname.split()[0]
        if readname not in mapped_reads:
            mapped_reads[readname]=[]
        cur_entry={}
        cur_entry['mapped_name']=mapped_name
        cur_entry['pos0']=pos0
        cur_entry['num_alignments']=num_alignments
        cur_entry['seq']=read
        cur_entry['strand']=strand
        mapped_reads[readname].append(cur_entry)
        readlist.add(readname)

    c=0
    cells=set()
    cell_dict={}
    with gzip.open(args.R1, "r") as handle:    
        for record in SeqIO.parse(handle, "fastq"):
            c+=1
            if c%1000000==0:
                print c
            cell_bc=str(record.seq[args.cell_bc_start:args.cell_bc_end])
            cells.add(cell_bc)
            if cell_bc not in cell_dict:
                cell_dict[cell_bc]=0
            cell_dict[cell_bc]+=1
            if record.id in mapped_reads:
                #cur_read_idx+=1
                #cur_read=readlist[cur_read_idx]
                cell_bc=str(record.seq[args.cell_bc_start:args.cell_bc_end])
                umi=str(record.seq[args.UMI_start:args.UMI_end])
                for entry in mapped_reads[record.id]:
                    out.write(entry['mapped_name']+'\t'+entry['seq']+'\t'+entry['pos0']+'\t'+entry['strand']+'\t'+cell_bc+'\t'+umi+'\t'+entry['num_alignments']+'\t'+record.id+'\n')

    print len(list(cells))
    for cell_bc in cell_dict:
        if int(cell_dict[cell_bc])>1000:
            out_cells.write(cell_bc+'\t'+str(cell_dict[cell_bc])+'\n')

main()
