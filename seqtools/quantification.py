class TPMCalculator:
   def __init__(self,txome):
      self.transcripts = {}
      for x in txome.transcripts:
         self.transcripts[x.name] = {}
         self.transcripts[x.name]['length'] = x.length
         self.transcripts[x.name]['count'] = 0
         self.transcripts[x.name]['TPM'] = 0
         self.transcripts[x.name]['FPKM'] = 0
         self.transcripts[x.name]['RPK'] = 0
      self._calculated = True
   def add_count(self,name,cnt=1):      
      self.transcripts[name]['count'] += cnt
      self._calculated = False
   def TPM(self,name):
      if not self._calculated: self.calculate()
      return self.transcripts[name]['TPM']
   def FPKM(self,name):
      if not self._calculated: self.calculate()
      return self.transcripts[name]['FPKM']
   def count(self,name):  return self.transcripts[name]['count']
   def calculate(self):
      """do the TPM calculation"""
      self._calculated = True
      for name in self.transcripts:
         self.transcripts[name]['RPK'] = (float(self.transcripts[name]['count'])/float(self.transcripts[name]['length']))/float(1000)
      tot = 0.0
      for name in self.transcripts:
         tot += self.transcripts[name]['RPK']
      tot = tot/float(1000000)
      for name in self.transcripts:
         self.transcripts[name]['TPM'] = self.transcripts[name]['RPK']/tot
      """do the FPKM calculation"""
      tot = 0
      for name in self.transcripts: tot += self.transcripts[name]['count']
      tot = float(tot)/float(1000000)
      for name in self.transcripts:
         if tot > 0:
            rpm = float(self.transcripts[name]['count'])/float(tot)
            self.transcripts[name]['FPKM'] = 1000*rpm/self.transcripts[name]['length']
