"""a module provide functions to conver transcripts to other formats

"""
import sys

def transcript_to_gpd_line(self,transcript_name=None,gene_name=None,strand=None):
    """Get the genpred format string representation of the mapping

    :param transcript_name:
    :param gene_name:
    :param strand:
    :type transcript_name: string
    :type gene_name: string
    :type strand: string
    :return: GPD line
    :rtype: string
    """
    self._initialize()
    tname = self._transcript_name
    gname = self._gene_name
    dir = self._direction
    # check for if we just have a single name
    if not tname and not gname:
      if self._name:
        tname = self._name
        gname = self._name
    if not tname: tname = transcript_name
    if not gname: gname = gene_name
    if not dir: dir = strand
    if not tname or not gname or strand:
      sys.stderr.write("ERROR:  transcript name and gene name and direction must be set to output a gpd line or use get_fake_gpd_line()\n")
    out = ''
    out += tname + "\t"
    out += gname + "\t"
    out += self.exons[0].rng.chr + "\t"
    out += dir + "\t"
    out += str(self.exons[0].rng.start-1) + "\t"
    out += str(self.exons[-1].rng.end) + "\t"
    out += str(self.exons[0].rng.start-1) + "\t"
    out += str(self.exons[-1].rng.end) + "\t"
    out += str(len(self.exons)) + "\t"
    out += str(','.join([str(x.rng.start-1) for x in self.exons]))+','+"\t"
    out += str(','.join([str(x.rng.end) for x in self.exons]))+','
    return out

def transcript_to_fake_psl_line(self,ref):
    """Convert a mapping to a fake PSL line

    :param ref: reference genome dictionary
    :type ref: dict()
    :return: psl line
    :rtype: string
    """
    self._initialize()
    e = self
    mylen = 0
    matches = 0
    qstartslist = []
    for exon in self.exons:
      mylen = exon.rng.length()
      matches += mylen
      qstartslist.append(matches-mylen)
    qstarts = ','.join([str(x) for x in qstartslist])+','
    oline =  str(matches)+"\t" # 1
    oline += "0\t" # 2
    oline += "0\t" # 3
    oline += "0\t" # 4
    oline += "0\t" # 5
    oline += "0\t" # 6
    oline += "0\t" # 7
    oline += "0\t" # 8
    oline += e.get_strand()+"\t" # 9
    oline += e.get_transcript_name()+"\t" # 10
    oline += str(matches)+"\t" # 11
    oline += "0\t" # 12
    oline += str(matches)+"\t" # 13
    oline += e.get_chrom()+"\t" # 14
    oline += str(len(ref[e.get_chrom()]))+"\t" # 15
    oline += str(e.exons[0].rng.start-1)+"\t" # 16
    oline += str(e.exons[-1].rng.end)+"\t" # 17
    oline += str(len(e.exons))+"\t" # 18
    oline += ','.join([str(e.exons[x].rng.end-(e.exons[x].rng.start-1)) for x in range(0,len(e.exons))])+','+"\t" # 19
    oline += qstarts + "\t" # 20
    oline += ','.join([str(x.rng.start-1) for x in e.exons])+',' # 21
    return oline
