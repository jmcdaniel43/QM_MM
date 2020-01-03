#  module load anaconda3/latest
#  source activate psi4_beta
#  export PSIDATADIR=~/.conda/envs/psi4_beta/share/psi4/

import psi4
import psi4.core as core
import subprocess

acnt=psi4.geometry("""
0  1
  C    0.000000    0.000000    0.000000
  C    0.056130    0.023070    1.447480
  N    0.101910    0.043860    2.594820
  H    0.235630   -0.997960   -0.368840
  H    0.720430    0.707160   -0.412730
  H   -0.998390    0.275760   -0.339710
symmetry c1
""")


#********* options for DMA ***************
# these options should already be in psi4_dma_datafile
#radii=['C',0.53,'N',0.53,'H',0.53]
#core.set_local_option('GDMA', 'GDMA_RADIUS', radii )
#core.set_local_option('GDMA', 'GDMA_LIMIT', 4 )
#core.set_local_option('GDMA', 'GDMA_MULTIPOLE_UNITS', 'AU' )
#core.set_local_option('GDMA', 'GDMA_SWITCH', 0.0 )

#*********  External Field
core.set_local_option('SCF','PERTURB_H',1)
core.set_local_option('SCF','PERTURB_WITH','DIPOLE')
core.set_local_option('SCF','PERTURB_DIPOLE',[ 0.0, 0.0, 0.0] )


( energy , wfn ) = psi4.energy('pbe0/def2-SVP', return_wfn=True )


# calling psi4.gdma will exit python, so we just repeat some of the code in this method to do it here
fw = core.FCHKWriter(wfn)
fchkfile = 'molecule' + '.fchk'
fw.write(fchkfile)

# do DMA, fit charges to DMA, read in charges
p = subprocess.Popen("./fit_charges.sh", stdout=subprocess.PIPE)
out, err = p.communicate()
charges=out.splitlines()
# for some reason these are byte encoded, convert back to string then float
charges=[ float(charges[j].decode("utf-8")) for j in range(len(charges)) ]

print( charges )





