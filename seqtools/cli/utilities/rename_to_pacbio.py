import sys,argparse
from seqtools.format.fasta import FASTAStream
from seqtools.format.fastq import FASTQStream

def main(args):
  input = sys.stdin
  if args.input=='-': input = sys.stdin
  else: input= open(args.input)
  output = sys.stdout
  if args.output: output = open(args.output,'w')
  if args.output_table:  args.output_table= open(args.output_table,'w')
  if args.gpd:
    z = 0
    for line in input:
        f = line.rstrip().split("\t")
        z+=1
        name = 'm150101_010101_11111_c111111111111111111_s1_p0/'+str(z)+'/ccs'
        if args.output_table: args.output_table.write(f[0]+"\t"+name+"\n")
        f[0] = name
        f[1] = name
        output.write("\t".join(f)+"\n")
    output.close()
    if args.output_table:
      args.output_table.close()
    return


  if args.fasta:
    input = FASTAStream(input)
  elif args.fastq:
    input = FASTQStream(input)
  z = 0
  for e in input:
    z+=1
    name = 'm150101_010101_11111_c111111111111111111_s1_p0/'+str(z)+'/ccs'
    if args.fastq:
      output.write( '@'+name+"\n"+ str(e.sequence)+"\n"+ '+'+"\n"+str(e.qual)+"\n")
    elif args.fasta:
      output.write('>'+name+"\n"+str(e.sequence)+"\n")
    if args.output_table: args.output_table.write(e.name+"\t"+name+"\n")
  output.close()
  if args.output_table: args.output_table.close()

def do_inputs():
  parser = argparse.ArgumentParser()
  parser.add_argument('input',help="Use - for STDIN")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--fasta',action='store_true')
  group.add_argument('--fastq',action='store_true')
  group.add_argument('--gpd',action='store_true')
  parser.add_argument('--output_table',help='save coversion to file')
  parser.add_argument('-o','--output')
  args = parser.parse_args()
  return args
def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
